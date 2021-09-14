set(0,'defaultFigureRenderer','painters')

% directory of preprocessed data. 
ppDir = '/mnt/mfs/selected_data/preprocessedEHS/preprocessedColumbiaUTAHARRAYSeizures'

% [20190612] starting with C7, since  I have those IEDs characterized. 
% looping over patients
nSzs = 14;
for ptsz = 6:nSzs
	close all; 
	clearvars -except ptsz ppDir nSzs 
	tic; fprintf('loading data...');
	switch ptsz
		case 1
			ptID = 'C5-1';
			load(fullfile(ppDir,'C5_MUAtimes-SHORT-1.mat'))
			load(fullfile(ppDir,'C5_1Kdownsampled_seizure-SHORT-1.mat'))
			
			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = 0; 
			Tend = 50; 
			timeWin = [Tstrt Tend];
			IWperiod = [125 140];
		case 2
			ptID = 'C5-2';
			load(fullfile(ppDir,'C5_MUAtimes-SHORT-2.mat'))
			load(fullfile(ppDir,'C5_1Kdownsampled_seizure-SHORT-2.mat'))

			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = 0; 
			Tend = 150; 
			timeWin = [Tstrt Tend];
			IWperiod = [125 150];
		case 3
			ptID = 'C5-3';
			load(fullfile(ppDir,'C5_MUAtimes-SHORT-3.mat'))
			load(fullfile(ppDir,'C5_1Kdownsampled_seizure-SHORT-3.mat'))

			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = 0; 
			Tend = 100; 
			timeWin = [Tstrt Tend];
			IWperiod = [120 160];
		case 4
			ptID = 'C7-1';
			% load MUA data 
			load(fullfile(ppDir,'C7_MUAtimes-SHORT.mat'))
			load(fullfile(ppDir,'C7_1Kdownsampled_seizure-SHORT.mat'))

			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = -12; 
			Tend = 12; 
			timeWin = [Tstrt Tend];
			IWperiod = [122 128];
		case 5
			ptID = 'CUCX2-1';
	 		% load MUA data 
			load(fullfile(ppDir,'CUCX2_sz1_MUA_times.mat'))
			load(fullfile(ppDir,'CUCX2_2Kdownsampled_seizure-1.mat'))

			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = -60; 
			Tend = 50; 
			timeWin = [Tstrt Tend];	
			IWperiod = [220 280];
		case 6
			ptID = 'CUCX2-2';
	 		% load MUA data 
			load(fullfile(ppDir,'CUCX2_sz2_MUA_times.mat'))
			load(fullfile(ppDir,'CUCX2_2Kdownsampled_seizure-2.mat'))

			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = 120; 
			Tend = 170; 
			timeWin = [Tstrt Tend];	
			IWperiod = [220 280]; 
		case 7
			ptID = 'CUCX3-1';
	 		% load MUA data 
			load(fullfile(ppDir,'CUCX3_sz1_MUA_times.mat'))
			load(fullfile(ppDir,'CUCX3_2Kdownsampled_seizure-1.mat'))

			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = 120; 
			Tend = 170; 
			timeWin = [Tstrt Tend];	
			IWperiod = [220 280];
		case 8
			ptID = 'CUCX3-2';
	 		% load MUA data 
			load(fullfile(ppDir,'CUCX3_sz1_MUA_times.mat'))
			load(fullfile(ppDir,'CUCX3_2Kdownsampled_seizure-1.mat'))

			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = 120; 
			Tend = 170; 
			timeWin = [Tstrt Tend];	
			IWperiod = [220 280];
		case 9
			ptID = 'CUCX3-3';
	 		% load MUA data 
			load(fullfile(ppDir,'CUCX3_sz3_MUA_times.mat'))
			load(fullfile(ppDir,'CUCX3_2Kdownsampled_seizure-3.mat'))

			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = 120; 
			Tend = 200; 
			timeWin = [Tstrt Tend];	
			IWperiod = [230 320];
		case 10
			ptID = 'CUCX3-4';
	 		% load MUA data 
			load(fullfile(ppDir,'CUCX3_sz4_MUA_times.mat'))
			load(fullfile(ppDir,'CUCX3_2Kdownsampled_seizure-4.mat'))

			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = 120; 
			Tend = 200; 
			timeWin = [Tstrt Tend];	
			IWperiod = [230 320];
		case 11 
			ptID = 'CUCX3-5';
	 		% load MUA data 
			load(fullfile(ppDir,'CUCX3_sz5_MUA_times.mat'))
			load(fullfile(ppDir,'CUCX3_2Kdownsampled_seizure-5.mat'))

			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = 120; 
			Tend = 200; 
			timeWin = [Tstrt Tend];	
			IWperiod = [230 320];
		case 12
			ptID = 'CUCX3-6';
	 		% load MUA data 
			load(fullfile(ppDir,'CUCX3_sz6_MUA_times.mat'))
			load(fullfile(ppDir,'CUCX3_2Kdownsampled_seizure-6.mat'))

			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = 120; 
			Tend = 200; 
			timeWin = [Tstrt Tend];	
			IWperiod = [220 310];
		case 13
			ptID = 'CUCX4-1';
	 		% load MUA data 
			load(fullfile(ppDir,'CUCX4_sz1_MUA_times.mat'))
			load(fullfile(ppDir,'CUCX4_2Kdownsampled_seizure-1.mat'))

			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = 120; 
			Tend = 200; 
			timeWin = [Tstrt Tend];	
			IWperiod = [220 300];
		case 14
			ptID = 'CUCX4-2';
	 		% load MUA data 
			load(fullfile(ppDir,'CUCX4_sz2_MUA_times.mat'))
			load(fullfile(ppDir,'CUCX4_2Kdownsampled_seizure-2.mat'))

			% pt specific details. 
			goodChannels = 1:96;
			Tstrt = 120; 
			Tend = 200; 
			timeWin = [Tstrt Tend];
			IWperiod = [220 310];
	end
	A = toc; fprintf('took %d seconds.\n',A)


	% generating firing rates with a slow kernel. 
	for ch = 1:length(goodChannels)
		MUAdata.times = mua_data.timestamps{goodChannels(ch)};
		MUAdataChs(goodChannels(ch)).times = mua_data.timestamps{goodChannels(ch)};

		% [20141211] PSTHs using chronux to determine direction of recruitment wave and potential direction of MUA firing per burst
		if isempty(psth(MUAdata,.250,'n',timeWin+120,1))
			channelPSTHs(goodChannels(ch),:) = zeros(1,length(channelPSTHs));
			channelPSTH10(goodChannels(ch),:) = zeros(1,length(channelPSTH10));
		else
			[channelPSTHs(goodChannels(ch),:),tPSTH,~] = psth(MUAdata,.5,'n',timeWin+120,1);
			[channelPSTH10(goodChannels(ch),:),tPSTH10,~] = psth(MUAdata,.010,'n',timeWin+120,1);
		end
	end 
	
	% zscoring to deal with large values. 
	% channelPSTHs = zscore(channelPSTHs,[],2);
	% [IWpeaks,IWtimeIdcs] = max(channelPSTHs(:,tPSTH>IWperiod(1) & tPSTH<IWperiod(2))');

	% determining timing ranks. 
	% [sortedIWtimes,IWranks] = sort(tPSTH(IWtimeIdcs));

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
	% RESPONSE = num2cell(tPSTH(IWtimeIdcs));
	% IWbetas = SpatialLinearRegression(RESPONSE,P,'switch_plot',0,'LossFun','L1');

	% accounting for different variable names. 
	if exist('seizure_downsampled','var')
		Seizure1K = seizure_downsampled; % it's really 2K. confirmed: line 21 of ~/Code/ /get_CUCX_seizures.m 
		clear seizure_downsampled
		dsFs = 2e3;
	else
		dsFs = 1e3;
	end

	% time windows. 
	if strcmp(ptID,'C7-1')
		plotWin = [115 130];
	else
		plotWin = [IWperiod(1)-5 IWperiod(2)+5];
	end
	
	% trimming data to fit into the plot win and upsampling. 
	tSec = linspace(0,(length(Seizure1K)/1000),length(Seizure1K));
	
%	LFP = Seizure1K(:,tSec>plotWin(1) & tSec<plotWin(2));
%	IW = resample(channelPSTHs(:,tPSTH>plotWin(1) & tPSTH<plotWin(2))',length(LFP),sum(tPSTH>plotWin(1) & tPSTH<plotWin(2)))';
%	BST = resample(channelPSTH10(:,tPSTH10>plotWin(1) & tPSTH10<plotWin(2))',length(LFP),sum(tPSTH10>plotWin(1) & tPSTH10<plotWin(2)))';
%	tSecWin = tSec(tSec>plotWin(1) & tSec<plotWin(2));

	% for looking at full sz tper.
	LFP = Seizure1K;
	IW = resample(channelPSTHs',length(LFP),length(tPSTH))';
	BST = channelPSTH10;
	tSecWin = tSec;

	clear Seizure1K channelPSTHs channelPSTH10 tSec

	% normalize firing rate between min and max
	IW = (IW-min(IW,[],2))./max(IW,[],2);
	BST = (BST-min(BST,[],2))./max(BST,[],2);

	% reshape matrices to tensors in array space. 
	for ix = 1:length(P)
		LFPmat(p1(ix),p2(ix),:) = LFP(ix,:); 
		IWmat(p1(ix),p2(ix),:) = IW(ix,:); 
		BSTmat(p1(ix),p2(ix),:) = BST(ix,:); 
	end
	LFPmat([1 1 10 10],[1 10 1 10],:) = NaN;
	IWmat([1 1 10 10],[1 10 1 10],:) = NaN;
	BSTmat([1 1 10 10],[1 10 1 10],:) = NaN;

	% TODO:: interpolate, but keep corners the same relative size, wihtout interpolation...
	% timing parameters. 
	nFrames = length(tSecWin);
	
	% colormap across channels for recruitment (black => early)
	cMap = colormap(hot(nFrames));

	% defining limits. 
	cLimsLFP = [min(min(LFP)) max(max(LFP))]
	cLimsAP = [0 max(max(BST))]
	cLimsIW = [0 max(max(IW))]

	% now making the movie. 
	movieDir = '/home/NIMASTER/ehs2149'
	movieFileName = fullfile(movieDir,sprintf('%s_ictalWavefrontMovie.avi',ptID));
	movieObj = VideoWriter(movieFileName,'Motion JPEG AVI')
	movieObj.FrameRate = 100;
	movieObj.Quality = 50;
	open(movieObj);

	% set up figure
	movH = figure(666);
	maximize(666);

	szAx = subplot(10,1,1:3);
	hold on
	plot(tSecWin,mean(LFP),'color','k')
	plot(tPSTH10,BST*100,'color','b');
	line([plotWin(1) plotWin(1)+3],[cLimsLFP(2)-300 cLimsLFP(2)-300],'color','k','linewidth',5)
	text(plotWin(1),cLimsLFP(2)+100,'3 seconds')
	hold off
	axis tight 

	% looping through time. 
	tic
	for f = 1   % :5:nFrames
		fprintf('plotted %d out of %d movie frames...(%.2f percent complete)',f,nFrames,100*(f/nFrames))

		% plot timing line.
		szAx = subplot(10,1,1:3);
		hold on
		line([tSecWin(f) tSecWin(f)],[cLimsLFP(2)-100 cLimsLFP(2)],'color','r','linewidth',2)
		hold off

		% plot action potential data across array. 
		subplot(10,3,[10 13 16 19 22 25 28]);
		imagesc(squeeze(LFPmat(:,:,f)),cLimsLFP)
		axis off square
		title('LFP')

		% plot IW data across the array.
		subplot(10,3,[11 14 17 20 23 26 29]);		
		imagesc(squeeze(IWmat(:,:,f)),cLimsIW)
		axis off square
		title('      unit firing: 250 ms smoothing (ictal wavefront)')

		% plot LFP data across array. 
		subplot(10,3,[12 15 18 21 24 27 30]);
		imagesc(squeeze(BSTmat(:,:,f)),cLimsAP)
		axis off square
		title('      unit firing: 10 ms smoothing (each burst)')

		colormap hot
	
		if isequal(f,1)
			saveas(666,sprintf('~/%s_videoPreview.pdf',ptID))
		end

		% saving data. 
		writeVideo(movieObj,getframe(movH))
	end
	A = toc;
	fprintf('plotting movie took %.2F minutes',A/60)
end
