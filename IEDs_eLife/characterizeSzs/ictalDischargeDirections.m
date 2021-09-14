set(0,'defaultFigureRenderer','painters')

% directory of preprocessed data. 
ppDir = '/mnt/mfs/selected_data/preprocessedEHS/preprocessedColumbiaUTAHARRAYSeizures'

% [20190612] starting with C7, since  I have those IEDs characterized. 
% looping over patients
nSzs = 14;
for ptsz = 1:nSzs

	clearvars -except ptsz ppDir nSzs 
	close all; 

	%%%%% all of the patient details %%%%%%
	tic; fprintf('loading data...');
	switch ptsz
		case 1
		end
	A = toc; fprintf('took %d seconds.\n',A)


	%%%%%%%%%%%%%% spatial map %%%%%%%%%%%%%%%%%
	if strcmp(ptID(1:4),'CUCX')
		% there is a different pinout for CUCX patients...
		if exist('goodChannels','var')
			[Ch,IA,IB] = intersect(electrodePinoutCUCX,goodChannels);
			[p1,p2] = ind2sub([10 10],IA);
			P = [p1,p2]*0.04; % convert to cm
		else
			[Ch,IA,IB] = intersect(electrodePinoutCUCX,1:96);
			[p1,p2] = ind2sub([10 10],IA);
			P = [p1,p2]*0.04; % convert to cm
		end
	else
		% ...than there was for the original C# patients.
		if exist('goodChannels','var')
			[Ch,IA,IB] = intersect(electrodepinout,goodChannels);
			[p1,p2] = ind2sub([10 10],IA);
			P = [p1,p2]*0.04; % convert to cm
		else
			[Ch,IA,IB] = intersect(electrodepinout,1:96);
			[p1,p2] = ind2sub([10 10],IA);
			P = [p1,p2]*0.04; % convert to cm
		end
	end	
	

	%%%%%%%%%%%% LFP voltage %%%%%%%%%%%%%%%%%
	% Loading the LFP, while accounting for different variable names. 
	if exist('seizure_downsampled','var')
		Seizure1K = seizure_downsampled; % it's really 2K. confirmed: line 21 of ~/Code/ /get_CUCX_seizures.m 
		clear seizure_downsampled
		dsFs = 2e3;
	else
		dsFs = 1e3;
	end
	tSec = linspace(0,length(Seizure1K)./dsFs,length(Seizure1K));

	% saving the LFP within the IW period and de-meaning
	LFP = Seizure1K(:,tSec>timeWin(1) & tSec<timeWin(2))-repmat(median(Seizure1K(:,tSec>timeWin(1) & tSec<timeWin(2)),2),1,sum(tSec>timeWin(1) & tSec<timeWin(2)));
	tSec = linspace(timeWin(1),timeWin(2),length(LFP));

	% filter in high frequency range. 
	[b,a] = butter(4,[80 200]./dsFs);
	for ch = 1:size(LFP,1)
		HG(ch,:) = filtfilt(b,a,LFP(ch,:));
	end

	%%%%%%%%%% Detect discharges %%%%%%%%%%%%%
	[bstAmps,bstTimeIdcs] = findpeaks(mean(smoothdata(abs(HG)','gaussian',25)'),'MinPeakProminence',10);

	% starting IDdata struct for saving data. 
	IDdata.burstTimesSamps = bstTimeIdcs;
	IDdata.burstTimesSecs = tSec(bstTimeIdcs);
	IDdata.burstAmplitudesHG = bstAmps;

	% loop over peaks
	nBsts = length(bstTimeIdcs)
	for bst = 2:nBsts
		% defining burst windows as 	
		if bst==1
			try
				binEdges = [tSec(bstTimeIdcs(bst)-min(diff(bstTimeIdcs))) tSec(bstTimeIdcs(bst+1)-min(diff(bstTimeIdcs)))];
			end
		elseif bst==nBsts
			try
				binEdges = [tSec(bstTimeIdcs(bst-1)+median(diff(bstTimeIdcs))) tSec(bstTimeIdcs(bst)+median(diff(bstTimeIdcs)))];
			end
		else
			try
				binEdges = [tSec(bstTimeIdcs(bst-1)+min(diff(bstTimeIdcs))) tSec(bstTimeIdcs(bst+1)-min(diff(bstTimeIdcs)))];
			end
		end

		% calculate ictal discharge travelin waves from voltage minima for comparison...
		if isequal(ptID,'C7-1')
			if (binEdges(2)-binEdges(1))<=0.001
				[negPeaks,peakLocs] = min(smoothdata(LFP(:,tSec>binEdges(1) & tSec<binEdges(1)+0.01)','gaussian',20));
			else
				[negPeaks,peakLocs] = min(smoothdata(LFP(:,tSec>binEdges(1) & tSec<binEdges(2))','gaussian',20));
			end
		else
			[negPeaks,peakLocs] = min(smoothdata(LFP(:,tSec>binEdges(1) & tSec<binEdges(2))','gaussian',20));
		end
			
		% spatial multilinear regression
		RESPONSE = num2cell(peakLocs/dsFs);
		[BSTbetas, ~, pVal] = SpatialLinearRegression(RESPONSE,P,'switch_plot',0,'LossFun','L1');
		V = pinv(BSTbetas(1:2));
	
		% that will be easy to incorporate into figure 2 scripts
		% so that it will be easy to visualize all direcional data on the same plot.
		IDdata.burstAmplitudesVpeak(bst) = nanmedian(negPeaks);
		IDdata.burstDirections(bst) = V(1) + V(2)*sqrt(-1);
		IDdata.burstSpeeds(bst) =  sqrt(sum(V.^2));
		IDdata.sigModel(bst) = pVal<0.05;

	end

	% removing speed outliers.
	speedOutliers = outliers(IDdata.burstSpeeds);
	IDdata.burstDirections(speedOutliers) = [];
	IDdata.burstSpeeds(speedOutliers) = [];
	IDdata.sigModel(speedOutliers) = [];
	IDdata.burstTimesSamps(speedOutliers) = [];
	IDdata.burstTimesSecs(speedOutliers) = [];
	IDdata.burstAmplitudesHG(speedOutliers) = [];
	IDdata.burstAmplitudesVpeak(speedOutliers) = [];

	if exist('preGenWin','var')
		IDdata.preGenWin = preGenWin; 
		IDdata.SecGenIdx = IDdata.burstTimesSecs>preGenWin(1) & IDdata.burstTimesSecs<preGenWin(2);
	end

	% colormap across channels for recruitment (black => early)
	cMapChs = colormap(copper(size(Seizure1K,1)));

	% visualize ictal discharges and seizures
	% in order to make sure that we're just detecting ictal discharges.
	burstFig = true;
	if burstFig
		% colormap across SIGNIFICANT bursts. ~ do open black circles for nonsignificant bursts and use this cmap for significant ones. 
		if ~isfield(IDdata,'SecGenIdx'); IDdata.SecGenIdx = true(1,length(IDdata.sigModel)); end
		cMapBst = colormap(cool(length(IDdata.burstDirections(IDdata.sigModel & IDdata.SecGenIdx))));

		%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		figure(1)
		
		% plot burst directions. 
		subplot(4,2,1)
		c = compass(IDdata.burstDirections(IDdata.sigModel & IDdata.SecGenIdx));
		for cl=1:length(c); c(cl).Color = cMapBst(cl,:); end
		
		% plot histograms of burst directions
		subplot(4,2,2)
		polarhistogram(angle(IDdata.burstDirections(IDdata.sigModel & IDdata.SecGenIdx)),18,'DisplayStyle','stairs')

		% plotting (mean?) voltage trace
		subplot(4,1,2)
		hold on
		plot(tSec,mean(LFP),'color','k')
		hold off
		axis tight
		xlim([timeWin])
		ylabel('')
		title(ptID)

		% plotting mean firing rate below voltage
		subplot(4,1,3)
		hold on
		plot(tSec,mean(smoothdata(abs(HG)','gaussian',50)'),'color','k')
		scatter(tSec(bstTimeIdcs(IDdata.sigModel & IDdata.SecGenIdx)),bstAmps(IDdata.sigModel & IDdata.SecGenIdx),2,cMapBst,'filled')
		scatter(tSec(bstTimeIdcs(~IDdata.sigModel)),bstAmps(~IDdata.sigModel),4,[0.3 0.7 0.3],'linewidth',0.5)
		hold off
		axis tight
		xlim([timeWin])
		ylim([0 8*std(mean(abs(HG)))])
		ylabel('high gamma power')
		xlabel('time (s)')

		% readout
		subplot(4,1,4)
		text(0,0,sprintf('%.2f percent significant travleing waves.',100*(sum(IDdata.sigModel & IDdata.SecGenIdx)/length(IDdata.sigModel & IDdata.SecGenIdx))))
		axis off

		% save figure
		halfMaximize(1,'left')
		saveas(1,fullfile('~/',sprintf('%s_IctalDischargeDirections.pdf',ptID)));

	end

	%%%%%%%%%%%%% SAVING DATA TO COMPARE WITH IEDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% saving ictal discharge data
	saveDir = '/storage/selected_data/elliotWorking/Data/IEDs/seizureData';
	save(fullfile(saveDir,sprintf('%s_IDdata.mat',ptID)),'IDdata','-v7.3')
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end







