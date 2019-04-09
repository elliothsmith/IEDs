% digging through IED detections to make a figure of each IED for visual inspection and get total numbers of IEDs across patients. 
parentDir = '/mnt/mfs/selected_data/elliotWorking/Data/IEDs';

set(0,'defaultfigurerenderer','painters')

%% [20190327] Notes: 
%	- This script is designed for Saturn2.columbia.neuro.edu
%	- This script should be run after IEDwaves_CU.m, where data about waves speeds should be saved.
%		I should now be able to load the wave data files for each IED and plot a figure for each IED. 

close all;

%% [20190327] dig through the IED detections for each patient
nPts = 8;
for pt = 3:nPts
    % (s)which patients
    clearvars -except parentDir pt nPts ALGO
    switch pt
        case 1
            ptID = 'c3';
            % sad magic variables resulting from my lack of prescience go here
            Fs = 3e4;
            negativeGoingIED = true;
        case 2
            ptID = 'c4';
            Fs = 3e4;
            negativeGoingIED = true;
        case 3
            ptID = 'c5';
            Fs = 3e4;
            negativeGoingIED = true;
			crashFiles = [2 19 34 47 48 49 53 54 56 64 71:73 156 157];
        case 4 
            ptID = 'c7';
            Fs = 3e4;
            negativeGoingIED = true;
			crashFiles = [2 50 51 96 118 132 139 140 147 148 153 154 155 180];
        case 5
            ptID = 'CUCX2';
            Fs = 3e4;
            negativeGoingIED = true;
    		crashFiles = [3 36 68 90 111 139 166 174 191] ;
    	case 6 
            ptID = 'CUCX3';
            Fs = 3e4;
            negativeGoingIED = true;
        	crashFiles = 2;
		case 7
            ptID = 'CUCX4';
            Fs = 3e4;
            negativeGoingIED = true;
			crashFiles = [2 11 19 29 38 39 44];
		case 8 
			ptID = 'CUCX5';
			Fs = 3e4;
			negativeGoingIED = true;
			crashFiles = 2; 
    end
    
    % where to save figs.
    figDir = '/mnt/mfs/selected_data/elliotWorking/Figs/IEDfigs/';
    IEDnumber = 0;

    % converting the 2D array map into a 96 X 2 matrix with cm units, P.
    map = electrodepinout;
    chs = sort(map(map>0));
    for i = numel(chs):-1:1
        [P(i,1),P(i,2)] = find(map == chs(i));
    end
    P = P*0.04; % cm
    
    % [20181017] looping over files in each IED directory
    dirList = dir(fullfile(parentDir,ptID));

	for fl = 1:length(dirList)
	fprintf('\npatient: %s, file number %d is %.2f megabytes',ptID,fl,dirList(fl).bytes./1e6)
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
                wavePointData = [IEDdata.fullBWdata(detectionIdcs(1),wavePointIdx).windowedData];
                
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

				%% This is where the other periodic error occurs. I should fix this...
                try
					if negativeGoingIED
                    	[localMinTimeIdcs,localMinAmps] = islocalmin(smoothedWPdata,2,'MinProminence',median(std(abs(smoothedWPdata),[],2)),'MaxNumExtrema',nChans*3);
                    	[absoluteMinAmps,absoluteMinTimeIdcs] = min(smoothedWPdata,[],2);
                	else
                    	[localMinTimeIdcs,localMinAmps] = islocalmax(smoothedWPdata,2,'MinProminence',median(std(abs(smoothedWPdata),[],2)),'MaxNumExtrema',nChans*3);
                    	[absoluteMinAmps,absoluteMinTimeIdcs] = max(smoothedWPdata,[],2);
                	end
					emptyFlag = false; 
                catch
					if ~isempty(smoothedWPdata)
						[localMinAmps,localMinTimeIdcs] = findpeaks(smoothedWPdata);
						[absoluteMinAmps,absoluteMinTimeIdcs] = min(smoothedWPdata,[],2);
						emptyFlag = false; 
					else
						fprintf('\ntraveling wave data is empty for this detection. \n skipping...');
						emptyFlag = true; 
					end
				end		

				if ~emptyFlag
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

                [IEDwaves(wp).beta, IEDwaves(wp).V, IEDwaves(wp).p] = ...
                            SpatialLinearRegression(RESPONSE,P,'switch_plot',0,'Lossfun','L1','Fs',Fs); %
                IEDwaves(wp).V = norm(pinv(IEDwaves(wp).beta(1:2))); %(in cm/s)

                %% visualize the data from this step
				figure(wp)
                  
			    % first plotting the full window.  	
				subplot(2,2,1)
				hold on
                plot(tSec,smoothedWPdata,'linewidth',0.5);
				scatter(negPeakTimes,negPeakValues,2,rgb('black'),'filled');
				hold off
                xlabel('time (s)')
                ylabel('LFP (uV)')
                xlim([0 1])
                axis square
                title(sprintf('smoothed data: %d ms gaussian kernel',(smoothWindowSize./Fs)*1000))
                    
				% sorting the colormap by timing
				[~,pIdcs] = sort(negPeakTimes,'ascend');
				tmpMap = colormap(copper(length(pIdcs)));
				cMap = tmpMap(pIdcs,:);

				% next plotting the zoomed-in waves. 
				subplot(2,2,3)
    			% colormap(copper);
		        zP = plot(tSec,smoothedWPdata,'linewidth',0.5);
				set(zP,{'Color'},num2cell(cMap,2));
				% replace colors
				% set(gca, 'ColorOrder', cMap, 'NextPlot', 'replacechildren');
				hold on
                scatter(negPeakTimes,negPeakValues,10,rgb('black'),'filled');
                hold off
                axis square tight
				xlim([median(negPeakTimes)-0.1 median(negPeakTimes)+0.1]);
                xlabel('time (s)')
                ylabel('LFP voltage')
                
				% plot the vectors
				subplot(2,2,2)
				c = compass(IEDwaves(wp).beta(1)+IEDwaves(wp).beta(2)*sqrt(-1));
				set(c,'Color',[0 0 0]),'linewidth',2
				xlabel('IED traveling wave direction')

				% plotted text
				subplot(2,2,4)
                title('IED traveling wave statistics:')
                hold on
                maxLim = 24;
                text(0,maxLim,sprintf('Subject: %s',ptID));
                text(0,maxLim-2,sprintf('model betas = [%d %d]',IEDwaves(wp).beta(1),IEDwaves(wp).beta(2)));
                if IEDwaves(wp).p<0.05
                    text(0,maxLim-4,sprintf('SIGNFICANT F-test against zero-slope: p=%.2d',IEDwaves(wp).p));
                else
                    text(0,maxLim-4,sprintf('NO DIFFERENCE against zero-slope: p=%.2d',IEDwaves(wp).p));
                end
                text(0,maxLim-6,sprintf('wave velocity = %.2d cm/s',IEDwaves(wp).V));
                hold off
                % deets
                ylim([0 25])
                axis off

				% saving figure
				saveas(wp,fullfile(figDir,sprintf('%s_IEDn%d_summary.pdf',ptID,IEDnumber)))
				close(wp)

				% save choice IED data
				if exist('/mnt/mfs/selected_data/elliotWorking/Data/IEDs/allIEDs.mat','file');
					load('/mnt/mfs/selected_data/elliotWorking/Data/IEDs/allIEDs.mat')
					sNum = length(thisIED)+1;
				else
					sNum = 1;
				end
				thisIED(sNum).patient = ptID;
				thisIED(sNum).file = fullfile(dirList(fl).folder,dirList(fl).name);
				thisIED(sNum).detectionSample = wavePoint;
				thisIED(sNum).minTimes = RESPONSE;
				thisIED(sNum).minAmps = localMinAmps;
				thisIED(sNum).waveBetas = IEDwaves(wp).beta;
				thisIED(sNum).waveVelo = IEDwaves(wp).V;
				thisIED(sNum).wavePval = IEDwaves(wp).p;
				save('/mnt/mfs/selected_data/elliotWorking/Data/IEDs/allIEDs.mat','thisIED','-v7.3')
				
				fprintf('\nfinished for IED %s',IEDnumber)
				end
			end % if statement for non detections
		end % for loop for each IED "wavepoint" [20190327] :: The fucking variable names in these scripts are not very informative. Makes it hard to come back to these analyses after leaving them for any period of time. This is only a problem since Cathy's and Tyler's servers went down, but good and intersting to know. 
	else
		fprintf('\nskipped %s',dirList(fl).name)

	end % if statement for skipping crashFiles
end % looping over files. 
end % main patient loop














%%%%%%% Some extra code that includes spectrograms if I have to go that route. 

				plotIEDwaves = false;
                if plotIEDwaves
                    % plotting results.
                    figure(wp)
                    plotmultipleaxes(1,2,2,0.1,wp);
                    imagesc(tSec,1:nChans,wavePointData);
                    axis square tight
                    xlabel('time (s)')
                    ylabel('UMA channels')
                    h = colorbar;
                    ylabel(h, 'LFP amplitude (uV)')
                    title(sprintf('IED # %d, all UMA channels',wp))
                    
                    plotmultipleaxes(3,2,2,0.1,wp);
                   
                    % [20181019] calculating IED spectrogram
                    for ch = 1:nChans
                        fPass = [1 200];
                        updateUser('calculating spectrograms',ch,20,nChans)
                        [W,period,~] = basewave4(wavePointData(ch,:),Fs,fPass(1),fPass(2),6,0);
                        Sft(:,:,ch) = abs((W))./repmat(1./logspace(0,log10(fPass(2)),length(period))',1,size(W,2));
                    end
                    
                    % determining spectral scales in Hz
                    scaleFreqs = 1./period;
                    
                    % plotting mean spectrogram across channels.
                    plotmultipleaxes(2,2,2,0.1,wp);
                    surf(tSec,scaleFreqs,squeeze(nanmean(Sft,3)),'edgecolor','none');
                    set(gca,'Yscale','log');
                    view(2);
                    axis square tight
                    xlabel('time (s)')
                    ylabel('frequency (Hz)')
                    h = colorbar;
                    ylabel(h, 'LFP power')
                    title('spectrogram across channels')
                    
                      
                    if ishandle(wp*1000)
                       % for the detection figure
                       halfMaximize(wp*1000,'left')
                       figName1000 = sprintf('%s/fig%d.pdf',tmpFigDir,wp*1000);
                       saveas(wp*1000,figName1000)
                    end
				end
				fprintf('\n finished for file number %d of %d (%s) (patient: %s) ',fl,length(dirList),dirList(fl).name,ptID)
 













