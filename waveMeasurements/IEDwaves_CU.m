% digging through IED detections to start quantifying traveling waves.
parentDir = '/mnt/mfs/selected_data/elliotWorking/Data/IEDs';

set(0,'defaultfigurerenderer','painters')

%% REMEMBER TO DELETE FIGURE PDFS!
% [20181120] built "deleteIEDPDFs.m" for this purpose.
%% OTHERWISE THE FIGURES CREATED HERE WILL BE ADDED TO THE PREVIOUS PDFS



%% TAKES A LONG TIME TO RUN ALL UTAH PATIENTS (~ 6 hours or so)
% 	This means that it will take days to weeks to run across the columbia patients. 
% 	started on 20181213::162000



%% which algorithm to use: 'mlinreg' or 'nonparam'
% ALGO = 'nonparam';



%% [20181017] to start: dig through the IED detections. (separate by patient?)
nPts = 8;
for pt = nPts
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
    tmpFigDir = '/mnt/mfs/selected_data/elliotWorking/Data/tmpIEDFigs/';
    
    % converting the 2D array map into a 96 X 2 matrix with cm units, P.
    map = electrodepinout;
    chs = sort(map(map>0));
    for i = numel(chs):-1:1
        [P(i,1),P(i,2)] = find(map == chs(i));
    end
    P = P*0.04; % cm
    
    % [20181017] looping over files in each IED directory
    dirList = dir(fullfile(parentDir,ptID));

	%% NOTE:: change the start for this loop variable if you have to run this again. 
	for fl = (crashFiles(end)+1):length(dirList)
	fprintf('\npatient: %s, file number %d is %.2f megabytes',ptID,fl,dirList(fl).bytes./1e6)
	if ~ismember(fl,crashFiles) & (dirList(fl).bytes./1e6)>200
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
        
        %% [20181018] looping over relevant detections and running
        % multilinear regression in order to determine IED traveling wave
        % speed and direction.
        nWavePoints = size(retainedDetectionIdcs,2);
        figList = {};
        for  wp = 1:nWavePoints
            wavePoint = median(allDetections(retainedDetectionIdcs(:,wp)));
            % loop to tell us which structure indices (channels) have the
            % corrcet data.
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
                
                % visualize the data from this step
                plotMins = true;
                if plotMins
                    figure(wp*100)
                    subplot(2,2,1)
                    plot(tSec,wavePointData)
                    xlabel('time (s)')
                    ylabel('LFP (uV)')
                    xlim([0 1])
                    axis square
                    title(sprintf('raw data: IED # %d',wp))
                    
                    subplot(2,2,2)
                    plot(tSec,smoothedWPdata)
                    xlabel('time (s)')
                    ylabel('LFP (uV)')
                    xlim([0 1])
                    axis square
                    title(sprintf('smoothed data: %d ms gaussian kernel',(smoothWindowSize./Fs)*1000))
                    
                    subplot(2,2,4)
                    hold on
                    plot(centers,weightingFunction.*N,'color',rgb('gray'))
                    scatter(tSecMat(localMinTimeIdcs),localMinAmps(localMinTimeIdcs)...
                        ,6,[0 0 0],'filled')
                    scatter(tSec(absoluteMinTimeIdcs),abs(absoluteMinAmps)...
                        ,4,rgb('lightseagreen'),'filled')
                    hold off
                    xlabel('time (s)')
                    ylabel('local minimum prominence')
                    xlim([0 1])
                    axis square
                    legend('weighted histogram','local minima','absolute minima')
                    title('minima detections: scatter')
                    % deets
                    xlabel('time (s)')
                    ylabel('peak prominence')
                    zlabel('local minima per bin')
                    axis xy
                    title('minima detections: histogram')
 %                   suptitle(sprintf('IED # %d',wp))
                end
                
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
                
                plotIEDwaves = true;
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
                    hold on
                    plot(tSec,smoothedWPdata);
                    scatter(negPeakTimes,negPeakValues,20,rgb('black'),'filled');
                    hold off
                    axis square tight
                    xlabel('time (s)')
                    ylabel('LFP voltage')
                    title(sprintf('smoothed data: %d ms gaussian kernel',(smoothWindowSize./Fs)*1000))
                    
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
                    
                    %% estimating travelling wave speed and direction.
%                     if strcmp(ALGO,'mlinreg')
                        % [20181017] multilinear regression step using least absolute deviation estimator
                        % why L1 regression? becuase it performed better re: MUA in Liou
                        % et al. (2017) JNeuEng
                        RESPONSE = num2cell(negPeakTimes');
                        figure(wp*1000)

                        [IEDwaves(wp).beta, IEDwaves(wp).V, IEDwaves(wp).p] = ...
                            SpatialLinearRegression(RESPONSE,P,'switch_plot',0,'Lossfun','L1','Fs',Fs); %
                        IEDwaves(wp).V = norm(pinv(IEDwaves(wp).beta(1:2))); %(in cm/s)
                        
                        if ishandle(wp*1000)
                            % for the detection figure
                            halfMaximize(wp*1000,'left')
                            figName1000 = sprintf('%s/fig%d.pdf',tmpFigDir,wp*1000);
                            saveas(wp*1000,figName1000)
                        end
                        
                        %[beta, V, p] = SpatialLinearRegression(T,P,'switch_plot',1) % Least square estimator
                        figure(wp)
                        plotmultipleaxes(4,2,2,0.2,wp)
                        % plotted text
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
                        
%                     elseif strcmp(ALGO,'nonparam')
                        % implementing nonparametric algorithm from Smith
                        % et al. (2016) NatComms
                        [~,channelOrder] = sort(negPeakTimes,'ascend');    % channelOrder: channels ordered by their burst minimae
                        
                        % burst orders as maps
                        Ich = zeros(length(channelOrder));
                        Jch = zeros(length(channelOrder));
                        
                        % looping over the bursts
                        for chz = 1:size(smoothedWPdata,1)
                            try
                                [I,J] = ind2sub(size(map),find(map==channelOrder(chz)));
                                % saving subscripts over Channels
                                Ich(chz) = I;
                                Jch(chz) = J;
                                
                                ArrayOrder(I,J) = chz;
                            catch
                                display(sprintf('couldn"t find channel %d in the map.',chz));
                            end                            
                        end
                        
                        %% calculating burst direction.
                        binaryImage = ArrayOrder~=0;
                        measurementsPlus = regionprops(logical(binaryImage),ArrayOrder,'WeightedCentroid');
                        measurementsMinus = regionprops(logical(binaryImage),size(smoothedWPdata,1)-ArrayOrder,'WeightedCentroid');
                        
                        Rminus = measurementsMinus.WeightedCentroid;
                        Rplus = measurementsPlus.WeightedCentroid;
                        
                        propDirection = Rplus-Rminus;
                        
                        [IEDwaves(wp).nonParamdirection,IEDwaves(wp).nonParamSpeed] = cart2pol(propDirection(1),propDirection(2));

%                     end
                end
				fprintf('\n finished for file number %d of %d (%s) (patient: %s) ',fl,length(dirList),dirList(fl).name,ptID)
            end
			end

         %% saving multi-page pdfs
        finalFigDir = sprintf('~/Figs/IEDs/%s/',ptID);
        if ~exist(finalFigDir,'dir')
            [~,MESSAGE,~] = mkdir(finalFigDir)
        end

         %    
         %    % handling figure stuff.
         %    if ishandle(wp) 
         %       % for the waves figure
         %        halfMaximize(wp,'right')
         %        figName = sprintf('%s/%s_wave%dProps',finalFigDir,dirList(fl).name(1:end-4),wp);
         %        print(wp,figName,'-dpdf')
		 %    	f = gcf
		 %    	f.Renderer
         %    end
         %    if ishandle(wp*100)
         %        % for the detection figure
         %        halfMaximize(wp*100,'left')
         %        figName100 =sprintf('%s/%s_wave%dDetection',finalFigDir,dirList(fl).name(1:end-4),wp);
         %        print(wp*100,figName100,'-dpdf')
		 %    	f = gcf
		 %    	f.Renderer
         %    end
			% 
           % if (ishandle(wp*100) && ishandle(wp) && exist('figName1000','var'))
           %     figList = cat(1,figList,{figName},{figName100},{figName1000});
           %     close(wp)
           %     close(wp*100)
                
           % elseif (ishandle(wp*100) && ishandle(wp))
           %     figList = {};
           % end

        end


       
		% Making a multi-page pdf. [Except this doesn't work]
%		finalFigName = [finalFigDir dirList(fl).name(1:end-4) '_allIEDs_travelingWaveSummary.pdf'];
%        append_pdfs(finalFigName, figList{:});
%        fprintf('\ncreated pdf of IED traveling waves: \n%s',finalFigName)
        
        %% saving udpated IED data structures with p values
        fprintf('\nsaving data...')
        tic
        save(fullfile(dirList(fl).folder,[dirList(fl).name(1:end-4) '_waves.mat']),'IEDwaves','-v7.3');
        A=toc;
        fprintf('\nsaving data took %d seconds...',A)
	else
		fprintf('\nskipped %s',dirList(fl).name)
	end % if statement for skipping crashFiles

    end % main file loop
end % main patient loop



%% NOTES::
% [20181017] NOTE:: There are several items I may have to add to
% may detection scripts if I want to realign the data or look at
% multiunits. The other option for doing this, is jsut going back
% to the original data, which may even be more efficient. Will have
% to thinkabout this a bit.



