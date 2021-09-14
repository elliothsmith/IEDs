% digging through IED detections to make a figure of each IED for visual inspection and get total numbers of IEDs across patients. 
parentDir = '/mnt/mfs/selected_data/elliotWorking/Data/IEDs';

set(0,'defaultfigurerenderer','painters')

%% [20190327] Notes: 
%	- This script is designed for Saturn2.columbia.neuro.edu
%	- This script should be run after IEDwaves.m, where data about waves speeds should be saved.

close all;

%% [20190327] dig through the IED detections for each patient
for pt = 1:7
    % (s)which patients
    clearvars -except parentDir pt ALGO
    switch pt
        case 1
            ptID = '';
            % sad magic variables resulting from my lack of prescience go here
            Fs = 3e4;
            negativeGoingIED = true;
			startFile = 1;
    end
	Fs = 3e4;

    % where to save figs.
    figDir = sprintf('/mnt/mfs/selected_data/elliotWorking/Figs/IEDfigs/%s/',ptID);  
	if ~exist(figDir,'dir'); mkdir(figDir); end

	% converting the 2D array map into a 96 X 2 matrix with cm units, P.
	map = electrodepinout;
	chs = sort(map(map>0));
	if strcmp(ptID,'c3'); chs = 1:79; end
	for i = numel(chs):-1:1
		[P(i,1),P(i,2)] = find(map == chs(i));
	end
	P = P*0.04; % cm
    
	if ~exist('startFile','var'); startFile = 1; end;
	if ~exist('IEDnumber','var'); IEDnumber = 0; end;
	if ~exist('crashFiles','var'); crashFiles = []; end;

	% [20181017] looping over files in each IED directory
    dirList = dir(fullfile(parentDir,ptID));

	for fl = startFile:length(dirList)
		fprintf('\npatient: %s, file number %d is %.2f megabytes',ptID,fl,dirList(fl).bytes./1e6)
		IEDSummaryFiles = dir(fullfile(parentDir,sprintf('allIEDs_%s_*.mat',ptID)));
		
		clearvars -except IEDnumber startFile crashFiles parentDir IEDSummaryFiles dirList P figDir negativeGoingIED ptID pt fl chs Fs
		
		if ~ismember(fl,crashFiles) & (dirList(fl).bytes./1e6)>200

			tic
			% loading preprocessed IED data.
			load(fullfile(dirList(fl).folder,dirList(fl).name));	

			%% [20181018] from the detection code. 
			detectionThreshold = 5;
			nSamps = IEDdata.parameters.downSamplingRate/2; % note this is double the previous window size (this is 1/2 second)
			nBins = ceil(IEDdata.resampledDataLength/nSamps);
			allDetections = sort(cat(1,IEDdata.detections.times));
			[detectionHisto,detectionEdges] = histcounts(allDetections,nBins);
			relevantLefts = detectionEdges([detectionHisto false]>detectionThreshold);
			relevantRights = detectionEdges([false detectionHisto]>detectionThreshold);
			
			% now find the detections within these bin limits.
			retainedDetectionIdcs = (allDetections>=relevantLefts & allDetections<=relevantRights);
			
	%		% loading IEDwave measurements. 
	%		load(fullfile(dirList(fl).folder,[dirList(fl).name(1:end-4) '_waves.mat']))
	%		A = toc; 
	%		fprintf('\nloading all IED data took %.2f seconds',A)

			%% [20181018] looping over relevant detections and running
			% multilinear regression in order to determine IED traveling wave
			% speed and direction.
			nWavePoints = size(retainedDetectionIdcs,2);
			for  wp = 1:nWavePoints
				IEDnumber = IEDnumber+1;
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
					try
						wavePointData = [IEDdata.fullBWdata(detectionIdcs(1),wavePointIdx).windowedData];
					catch
						wavePointIdx = floor(length(wavePointIdx./2));
						wavePointData = [IEDdata.fullBWdata(detectionIdcs(1),wavePointIdx).windowedData];
					end

					% extracting multiunit activity
					if ~isempty(wavePointData)
						[mua,muaTimes,Threshold] = IED_MUA(wavePointData);
					else
						 muaTimes = cell(1,length(chs));
					end

					% detect maximal descent
					% [20181024] very few sgnificant waves with this method
					[maxDs,maxDtimings] = min(diff(wavePointData,1,2),[],2);
					
					% params
					nChans = size(wavePointData,1);
					% [20181026] in case there are more channels in the data
					% than in the electrode map
					if gt(nChans,numel(chs))
						nChans = numel(chs);
						wavePointData = wavePointData(1:nChans,:);
					end
					tSec = linspace(0,IEDdata.parameters.windowSize,size(wavePointData,2));
					tSecMat = repmat(tSec,nChans,1);
				 
					% [20181024] now looking at many local minima and plotting
					% points used in multilinear regression model.
					binFactor = 40;
					smoothWindowSize = Fs./binFactor;
					smoothedWPdata = smoothdata(wavePointData,2,'gaussian',smoothWindowSize);

					if negativeGoingIED
						[localMinTimeIdcs,localMinAmps] = islocalmin(smoothedWPdata,2,'MinProminence',median(std(abs(smoothedWPdata),[],2)),'MaxNumExtrema',nChans*3);
						[absoluteMinAmps,absoluteMinTimeIdcs] = min(smoothedWPdata,[],2);
					else
						[localMinTimeIdcs,localMinAmps] = islocalmax(smoothedWPdata,2,'MinProminence',median(std(abs(smoothedWPdata),[],2)),'MaxNumExtrema',nChans*3);
						[absoluteMinAmps,absoluteMinTimeIdcs] = max(smoothedWPdata,[],2);
					end


						% minimum detection parameters.
						edges = 0:0.02:1;
						weightingFunction = [0.1*ones(1,(floor(length(edges)/2))-1) linspace(1,0.5,ceil(length(edges)/2))];
						[N,edges,bindices] = histcounts(tSecMat(localMinTimeIdcs),edges);
						centers = (edges(1:end-1) + edges(2:end))/2;

						%% [20181024] finding the local minima in the largest
						% weighted histogram bin
						[~,maxHisto] = max(weightingFunction.*N);
						% which chans have local minima in the range of interest
						if isequal(maxHisto,1)
							maxHisto = 25;
						elseif isequal(maxHisto,50)
							maxHisto = 25;
						end
						minChans = logical(sum(localMinTimeIdcs(:,tSec>centers(maxHisto-1) & tSec<centers(maxHisto+1)),2));
						bindicesStartPoint = find(diff(tSec<centers(maxHisto-1)));
						% calculating mins over all channels wihtin the bin.
						[minsAllChans,minTiming] = min(smoothedWPdata(:,tSec>centers(maxHisto-1) & tSec<centers(maxHisto+1)),[],2);
					
						%% [20181024] these are the important variables.~~~~~~~~~~~
						negPeakTimes = tSec(bindicesStartPoint+minTiming);
						negPeakValues = min(smoothedWPdata(:,bindicesStartPoint+minTiming),[],2);
						%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					
						% [20181017] multilinear regression step using least absolute deviation estimator
						% why L1 regression? becuase it performed better re: MUA in Liou
						% et al. (2017) JNeuEng
						RESPONSE = num2cell(negPeakTimes');

						% linear regression for LFP
						[IEDwaves(wp).beta, IEDwaves(wp).V, IEDwaves(wp).p] = ...
								SpatialLinearRegression(RESPONSE,P,'switch_plot',0,'Lossfun','L1','Fs',Fs); %
						%  IEDwaves(wp).V = norm(pinv(IEDwaves(wp).beta(1:2))); %(in cm/s)

						% linear regression for spikes
						try
							[IEDwaves(wp).spike_beta, IEDwaves(wp).spike_V, IEDwaves(wp).spike_p] = ...
								SpatialLinearRegression(muaTimes,P,'switch_plot',0,'Lossfun','L1','Fs',Fs); %
						catch % in case the model is poorly determined. 
							IEDwaves(wp).spike_beta = [NaN; NaN; NaN]; 
							IEDwaves(wp).spike_V = NaN; 
						   	IEDwaves(wp).spike_p = NaN; 
						end

                        fNum = length(IEDSummaryFiles);
                        if exist(sprintf('/mnt/mfs/selected_data/elliotWorking/Data/IEDs/allIEDs_%s_%d.mat',ptID,fNum),'file')
                         % load(sprintf('/mnt/mfs/selected_data/elliotWorking/Data/IEDs/allIEDs_%s_%d.mat',ptID,fNum))
	                        sNum = length(thisIED)+1;
	                    else
	                        sNum = 1;
	                    end

	                     thisIED(sNum).patient = ptID;
                        thisIED(sNum).file = fullfile(dirList(fl).folder,dirList(fl).name);
                        thisIED(sNum).detectionSample = wavePoint;
                        thisIED(sNum).waveBetas = IEDwaves(wp).beta;
                       thisIED(sNum).waveVelo = IEDwaves(wp).V;
                        thisIED(sNum).wavePval = IEDwaves(wp).p;
                         thisIED(sNum).waveBetasMUA = IEDwaves(wp).spike_beta;
                       thisIED(sNum).waveVeloMUA = IEDwaves(wp).spike_V;
                        thisIED(sNum).wavePvalMUA = IEDwaves(wp).spike_p;
                       thisIED(sNum).IEDnum = IEDnumber;
                        thisIED(sNum).medianWF = median(wavePointData);
                        thisIED(sNum).muaTimes = muaTimes;
                        thisIED(sNum).muaThresh = Threshold;
						thisIED(sNum).medAmp = median(negPeakValues);
						thisIED(sNum).meanAmp = mean(negPeakValues);
                        save(sprintf('/mnt/mfs/selected_data/elliotWorking/Data/IEDs/allIEDs_%s_%d.mat',ptID,fNum),'thisIED','-v7.3')

                        fprintf('\nfinished for IED %d (file number %d).',IEDnumber,fl)


					end
			end % for loop for each IED "wavepoint" [20190327] :: The fucking variable names in these scripts are not very informative. Makes it hard to come back to these analyses after leaving them for any period of time. This is only a problem since Cathy's and Tyler's servers went down, but good and intersting to know. 
		else
		fprintf('\nskipped %s',dirList(fl).name)

		end % if statement for skipping crashFiles
		clearvars -global
	end % looping over files. 
end % main patient loop

quit;








