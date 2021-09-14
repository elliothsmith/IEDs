

set(0,'defaultFigureRenderer','painters')

% directory of preprocessed data. 
ppDir = '/mnt/mfs/selected_data/preprocessedEHS/preprocessedColumbiaUTAHARRAYSeizures'
% looping over patients
nSzs = 14;
for ptsz = [1:6 13:14]

	clearvars -except ptsz ppDir nSzs 

	tic; fprintf('loading data...');
	switch ptsz
		case 1
	end
	A = toc; fprintf('took %d seconds.\n',A)

	% generating firing rates with a slow kernel. 
	for ch = 1:length(goodChannels)
		MUAdata.times = mua_data.timestamps{goodChannels(ch)};
		MUAdataChs(goodChannels(ch)).times = mua_data.timestamps{goodChannels(ch)};

		% [20141211] PSTHs using chronux to determine direction of recruitment wave and potential direction of MUA firing per burst
		if isempty(psth(MUAdata,.5,'n',timeWin+120,1))
			channelPSTHs(goodChannels(ch),:) = zeros(1,length(channelPSTHs));
			channelPSTH10(goodChannels(ch),:) = zeros(1,length(channelPSTH10));
		else
			[channelPSTHs(goodChannels(ch),:),tPSTH,~] = psth(MUAdata,.5,'n',timeWin+120,1);
			[channelPSTH10(goodChannels(ch),:),tPSTH10,~] = psth(MUAdata,.010,'n',timeWin+120,1);
		end
	end 
	
	% zscoring to deal with large values. 
	channelPSTHs = zscore(channelPSTHs,[],2);
	[IWpeaks,IWtimeIdcs] = max(channelPSTHs(:,tPSTH>IWperiod(1) & tPSTH<IWperiod(2))');

	% determining timing ranks. 
	[sortedIWtimes,IWranks] = sort(tPSTH(IWtimeIdcs));

	% spatial map. 
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
	
	% spatial multilinear regression
	plotRegression = true;
	if plotRegression
		RESPONSE = num2cell(tPSTH(IWtimeIdcs));
		[IWbetas,~,pVal] = SpatialLinearRegression(RESPONSE,P,'switch_plot',1,'LossFun','L1');
		saveas(gcf,sprintf('~/%s_IWregression.pdf',ptID))
	else
		RESPONSE = num2cell(tPSTH(IWtimeIdcs));
		[IWbetas,~,pVal] = SpatialLinearRegression(RESPONSE,P,'switch_plot',0,'LossFun','L1');
	end
	% colormap across channels for recruitment (black => early)
	cMap = colormap(copper(length(IWranks)));

	% accounting for different variable names. 
	if exist('seizure_downsampled','var')
		Seizure1K = seizure_downsampled; % it's really 2K. confirmed: line 21 of ~/Code/ /get_CUCX_seizures.m 
		clear seizure_downsampled
		dsFs = 2e3;
	else
		dsFs = 1e3;
	end

	% convert RESPONSE to matrix
	for ix = 1:length(RESPONSE); respMat(p1(ix),p2(ix)) = RESPONSE{ix}; rankMat(p1(ix),p2(ix)) = IWranks(ix);  end; 
	respMat([1 1 10 10],[1 10 1 10]) = NaN;

	% time windows. 
	if strcmp(ptID,'C7-1')
		plotWin = [115 130];
	else
		plotWin = [115 145];
	end
	plotWin = [IWperiod(1)-40 IWperiod(2)+50];

	plotFigure3 = true; 
	if plotFigure3
		%% FIGURE 3 STARTS HERE !!!
		% plotting a raster plot
		figure(1)
		% plotting (mean?) voltage trace
		subplot(6,2,1)
		plot(linspace(0,length(Seizure1K)./dsFs,length(Seizure1K)),mean(Seizure1K),'color','k')
		axis tight off
		xlim(plotWin)	
		title(ptID)

		% plotting mean firing rate below voltage
		subplot(6,2,3)
		plot(tPSTH10,mean(channelPSTH10),'color',[0.7 0.7 0.7])
		axis tight off
		xlim(plotWin)

		% TODO:: figure out what to put in this panel. 
		subplot(3,2,2)
		text(0,5,sprintf('IW speed: %.2f',abs(pinv(IWbetas(1:2)))))
		text(0,3,sprintf('model p-value: %.2f',pVal))

		% plotting raster
 		subplot(3,2,3)
 		plotSpikeRaster(mua_data.timestamps(IWranks)','VertSpikeHeight',2);
 		% deets
 		axis tight square
 		xlim(plotWin)
 		xlabel('time (s)')
 		ylabel('channels')

		% plotting firing rate signatures of ictal recruitment
		subplot(3,2,5)
		cp = plot(tPSTH,channelPSTHs(IWranks,:));
		set(cp,{'Color'},num2cell(cMap,2));
		% deets
		axis tight square
		xlim(plotWin)
		xlabel('time (s)')
		ylabel('per-channel firing rates (spikes/s)')
		
		% plotting regression plane
		subplot(3,2,4)
		hold on
		imagesc(respMat)
	% 	scatter3(p1,p2,cell2mat(RESPONSE),10,[0 0 0],'filled');
	% 	[P1U,P2U] = meshgrid(sort(unique(p1)),sort(unique(p2)));
	%	f = scatteredInterpolant(p1,p2,P*IWbetas(1:2) + IWbetas(3));
	%	Z = f(P1U,P2U);
	%	mesh(P1U,P2U,Z,'facealpha',0);
		hold off
		colormap(cMap)
		axis xy square tight off
		xlabel('microelectrodes')
		ylabel('microelectrodes')
		h = colorbar;
		ylabel(h,'ictal wavefront timing (s)')

		% plotting compass 
		subplot(3,2,6)
		V = pinv(IWbetas(1:2));
		cc = compass(V(1) + sqrt(-1)*V(2));
		title('direction of ictal wavefront')
		xlabel('space relative to array')

		% saving statistics
		IWdata.R = RESPONSE;
		IWdata.betas = IWbetas;
		IWdata.V = V;
		IWdata.pVal = pVal; 
		save(sprintf('/storage/selected_data/elliotWorking/Data/IEDs/seizureData/%s_IWbetas.mat',ptID),'IWdata','-v7.3')

		% saving figure
		halfMaximize(1,'left')
		saveDir = sprintf('~/Figs/IEDs/');
		saveas(1,fullfile(saveDir,sprintf('%s_szPlot.pdf',ptID)))
	end
end
