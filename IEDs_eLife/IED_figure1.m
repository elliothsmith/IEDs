% makes figure 1

%%%%%%%%%%% PICK AN IED %%%%%%%%%%%
% [20190613] at this point, the script will plot features of IEDs in the same file. 
IEDnumber = randperm(length(thisIED),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loading data for particular IED. 
load(thisIED(IEDnumber).file)

A = toc; 
fprintf('\nFinished loading data in %.2f seconds.',A)

% spatial map.
if isequal(length(patientID),2)
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
elseif strcmp(patientID(1:4),'CUCX')
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
end

%% [20181018] from the detection code.
detectionThreshold = 8;
nSamps = IEDdata.parameters.downSamplingRate/2; % note this is double the previous window size (this is 1/2 second)
nBins = ceil(IEDdata.resampledDataLength/nSamps);
allDetections = sort(cat(1,IEDdata.detections.times));
[detectionHisto,detectionEdges] = histcounts(allDetections,nBins);
relevantLefts = detectionEdges([detectionHisto false]>detectionThreshold);
relevantRights = detectionEdges([false detectionHisto]>detectionThreshold);

% now find the detections within these bin limits.
retainedDetectionIdcs = (allDetections>=relevantLefts & allDetections<=relevantRights);

% the wave points. 
wavePoints = size(retainedDetectionIdcs,2);

%% [20181018] looping over relevant detections and running
% multilinear regression in order to determine IED traveling wave
% speed and direction.
% [20190702] I guess I have to do this because the IED number info doesn't specify which detection...
nWavePoints = size(retainedDetectionIdcs,2);
for  wp = 1:nWavePoints
	IEDnumber = IEDnumber+1; % [20190702] why did I have this? it appears to just be for saving. 
	wavePoint = median(allDetections(retainedDetectionIdcs(:,wp)));
	
	% loop to tell us which structure indices (channels) have the
	% correct data.
	for dt = 1:length(IEDdata.detections)
		% which CHANNELS inlcude the chosen wavepoint
		dataDetections(dt) = ismember(wavePoint,IEDdata.detections(dt).times);
	end
	
	% [20181019] find the wavepoint index in the detections for the chosen channel.
	detectionIdcs = find(dataDetections);
	if isempty(detectionIdcs)
		figure(wp)
		tH = text(0,0,'NO DETECTIONS?');
		tH.FontWeight = 'bold';
		tH.FontSize = 24;
	else
		[~,wavePointIdx] = min(abs(IEDdata.detections(detectionIdcs(1)).times-wavePoint));
	
		% [20181019] This is all of the full-bandwidth data across channels for the particular IED.
		wavePointData = [IEDdata.fullBWdata(detectionIdcs(1),wavePointIdx).windowedData];
		Fs = IEDdata.parameters.OGsamplingRate;

		% detection multiunit acativity. 
		[mua,muaTimes,Threshold] = IED_MUA(wavePointData);

		% params
		nChans = size(wavePointData,1);

		tSec = linspace(0,IEDdata.parameters.windowSize,size(wavePointData,2));
		tSecMat = repmat(tSec,nChans,1);

		% [20181024] now looking at many local minima and plotting
		% points used in multilinear regression model.
		binFactor = 40;
		smoothWindowSize = Fs./binFactor;
		smoothedWPdata = smoothdata(wavePointData,2,'gaussian',smoothWindowSize);

		% finding mins...
		[localMinTimeIdcs,localMinAmps] = islocalmin(smoothedWPdata,2,'MinProminence',median(std(abs(smoothedWPdata),[],2)),'MaxNumExtrema',nChans*3);
		[absoluteMinAmps,absoluteMinTimeIdcs] = min(smoothedWPdata,[],2);

		% minimum detection parameters.
		edges = 0:0.02:1;
		weightingFunction = [0.1*ones(1,(floor(length(edges)/2))-1) linspace(1,0.5,ceil(length(edges)/2))];
		[N,edges,bindices] = histcounts(tSecMat(localMinTimeIdcs),edges);
		centers = (edges(1:end-1) + edges(2:end))/2;
		
		%% [20181024] finding the local minima in the largest weighted histogram bin
		[~,maxHisto] = max(weightingFunction.*N);
		% which chans have local minima in the range of interest
		if isequal(maxHisto,1);	maxHisto = 25;	elseif isequal(maxHisto,50); maxHisto = 25;	end;
		minChans = logical(sum(localMinTimeIdcs(:,tSec>centers(maxHisto-1) & tSec<centers(maxHisto+1)),2));
		bindicesStartPoint = find(diff(tSec<centers(maxHisto-1)));
		% calculating mins over all channels wihtin the bin.
		[minsAllChans,minTiming] = min(smoothedWPdata(:,tSec>centers(maxHisto-1) & tSec<centers(maxHisto+1)),[],2);
															
		%% [20181024] these are the important variables.~~~~~~~~~~~
		negPeakTimes = tSec(bindicesStartPoint+minTiming);
		negPeakValues = min(smoothedWPdata(:,bindicesStartPoint+minTiming),[],2);
		%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		%% part to calculate travleing wave speeds. 
		RESPONSE = num2cell(negPeakTimes');
		[IEDwaves(wp).beta, IEDwaves(wp).V, IEDwaves(wp).p] = SpatialLinearRegression(RESPONSE,P,'switch_plot',0,'Lossfun','L1','Fs',Fs); %  IEDwaves(wp).V = norm(pinv(IEDwaves(wp).beta(1:2))); %(in cm/s)

		%% plot figure 1
		figure(wp)
        borderSize = 0.05;

		% sorting the colormap by timing
		[~,pIdcs] = sort(negPeakTimes,'ascend');
		cMap = colormap(copper(length(pIdcs)));		% first plotting the full window.

		% plot example discharage. 
		ax1 = plotmultipleaxes(1,2,3,borderSize,wp);
		hold on
		plot(tSec,smoothedWPdata,'linewidth',0.5,'color',[0.5 0.5 0.5]);
		plot(tSec,mean(smoothedWPdata),'linewidth',1,'color',[1 1 1]);
		% scatter(negPeakTimes,negPeakValues,2,rgb('black'),'filled');
		hold off
		xlabel('time (s)')
		ylabel('LFP (uV)')
		xlim([0 1])
		axis square
	
		% next plotting the zoomed-in waves.
		ax2 = plotmultipleaxes(4,2,3,borderSize,wp);
		zP = plot(tSec,smoothedWPdata,'linewidth',0.5);
		set(zP,{'Color'},num2cell(cMap,2));
		% replace colors
		% set(gca, 'ColorOrder', cMap, 'NextPlot', 'replacechildren');
		hold on
		scatter(negPeakTimes,negPeakValues,10,rgb('black'),'filled');
		hold off
		axis square tight
		xlim([median(negPeakTimes)-0.05 median(negPeakTimes)+0.05]);
		xlabel('time (s)')
		ylabel('LFP voltage')

		% [20181019] calculating IED spectrogram
		for ch = 1:nChans
			fPass = [1 250];
			updateUser('calculating spectrograms',ch,20,nChans)
			[W,period,~] = basewave4(wavePointData(ch,:),Fs,fPass(1),fPass(2),6,0);
			Sft(:,:,ch) = abs((W))./repmat(1./logspace(0,log10(fPass(2)),length(period))',1,size(W,2));
		end
		
		% determining spectral scales in Hz
		scaleFreqs = 1./period;
		
		% plotting mean spectrogram across channels.
		ax3 = plotmultipleaxes(3,2,3,borderSize,wp);
		surf(tSec(1:30:end),scaleFreqs,squeeze(nanmean(Sft(:,1:30:end,:),3)),'edgecolor','none');
		set(gca,'Yscale','log');
		view(2);
		axis square tight
		colormap(ax3,turbo)
		xlabel('time (s)')
		ylabel('frequency (Hz)')
		h = colorbar;
		ylabel(h, 'LFP power')
		
		% plot firing rates across the array (colored by time). 
		plotmultipleaxes(2,2,3,borderSize,wp)
		plotSpikeRaster(muaTimes,'PlotType','vertline');
		axis square tight
		xlim([median(negPeakTimes)-0.05 median(negPeakTimes)+0.05]);
		xlabel('time (samples)')
		ylabel('UMA channels')

		% plot a spatial depiction of the wave:
		for ix = 1:length(RESPONSE); respMat(p1(ix),p2(ix)) = RESPONSE{ix}; end;
		% making the corners NaNs. 
		respMat(1,1) = NaN; respMat(1,10) = NaN; respMat(10,1) = NaN; respMat(10,10) = NaN;

		%  plotting 
		ax4 = plotmultipleaxes(5,2,3,borderSize,wp);
		imagesc(respMat);
		colormap(ax4,copper)
		axis xy square off
		xlabel('microelectrodes')
		ylabel('microelectrodes')
		h = colorbar;
		ylabel(h,'IED minima timing (s)')
		
		% plot the wave vectors
		plotmultipleaxes(6,2,3,0.03,wp)
		V = pinv(IEDwaves(wp).beta(1:2));
		c = compass(V(1)+V(2)*sqrt(-1));
		set(c,'Color',[0 0 0],'linewidth',1)
		xlabel('IED traveling wave direction')

		fprintf('\nfinished plotting. now saving.')

		% saving figure
		tic
		saveas(wp,sprintf('~/%s_file%d_IEDn%d_Figure1_summary.pdf',patientID,fl,IEDnumber))
		close(wp)
		A = toc;
		fprintf('\nsaving figure took %d seconds',A)
	end % if statement checking for empty deteection data
end % looping over retained detections. 





