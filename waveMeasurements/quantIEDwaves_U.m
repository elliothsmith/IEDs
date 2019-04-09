% digging through IED detections to start quantifying traveling waves.
%% This script is fast if not calculating waves. Slower if recaluclating vars.
%% [20181206] make sure you've used the distribution version of Jyunyou's
%% multilinear regression code to calculate traveling wave betas and velocity.

% defining the directory where patientID subdirectories are located.
parentDir = '~/data/IEDs';


%% [20181017] to start: dig through the IED detections. (separate by patient?)
nPts = 5;
IEDwaveCount = zeros(1,nPts); % number of detected IEDs overall.
for pt = [1 3] % 1:nPts
    % (s)which patients
    clearvars -except parentDir pt nPts sigIEDwaveVelocity sigIEDwaveDirection
    switch pt
        case 1
            ptID = 'nu';
            % sad magic variables resulting from my lack of prescience go here
            Fs = 3e4;
            ptHasSz = true;
            % seizure files and times.
            seizureFiles = {'NU_MUAtimes-SHORT-sz1.mat','NU_MUAtimes-SHORT-sz2.mat'};
            preictalTimes = [117 121];
            seizureDurations = [43 31]; % in s.
            IWwins = [0 9;0 15]; % in tPSTH time
            goodChannels = [1 3:11 13:37 39:64 68:69 73:96];
            deleteWaveFiles = false;
        case 2
            ptID = 'gamma';
            Fs = 3e4;
            ptHasSz = false;
            % seizure files and times.
            goodChannels = 1:96;
            deleteWaveFiles = false;
        case 3
            ptID = 'epsilon';
            Fs = 3e4;
            ptHasSz = true;
            seizureFiles = {'epsilon_MUAtimes-1.mat','epsilon_MUAtimes-2.mat'};
            preictalTimes = [121 121];
            seizureDurations = [75 315]; % in s.
            IWwins = [0 12; 60 90]; % in tPSTH time            
            
            goodChannels = 1:96;
            deleteWaveFiles = false;
        case 4
            ptID = 'eta';
            Fs = 3e4;
            ptHasSz = false;
            goodChannels = 1:96;
            deleteWaveFiles = true;
        case 5
            ptID = 'zeta';
            Fs = 3e4;
            ptHasSz = false;
            goodChannels = 1:96;
            deleteWaveFiles = true;
    end
    %% for deleting or saving calulations for all patients.
%     deleteWaveFiles = true;
    
    
    % where to save figs.
    tmpFigDir = '/media/user1/data4TB/data/IEDs/tmpIEDFigs';
    
    % converting the 2D array map into a 96 X 2 matrix with cm units, P.
    map = electrodepinout; % assuming that all of the utah utah arrays had the same pinout, which is the same as those from the early columbia cases.
    
    % first way of calculating P
    if exist('goodChannels','var')
        [Ch,IA,IB] = intersect(map,goodChannels);
        [p1,p2] = ind2sub([10 10],IA);
        P = [p1,p2]*0.04; % cm
    else
        [Ch,IA,IB] = intersect(map,1:96);
        [p1,p2] = ind2sub([10 10],IA);
        P = [p1,p2]*0.04; % cm
    end
    
    %     % second way of calculating P. [20181205] they correspond!
    %     ch = sort(map(map>0));
    %     for i = numel(ch):-1:1
    %         [P(i,1),P(i,2)] = find(map == ch(i));
    %     end
    %     P2 = P*0.04; % cm
    
    % pre-saving data to save time.
    ptWaveFile = sprintf('/media/user1/data4TB/data/IEDs/%s/%s_waveSpeedsAndDirections.mat',ptID,ptID);
    
    
    %% in case you want to start over with these analyses.
    % [20181204] takes about 20 mins to run all of them.
    if (deleteWaveFiles && exist(ptWaveFile,'file'))
        delete(ptWaveFile)
    end
    
    
    % accessing saved, preprocessed IED data.
    dirList = dir(fullfile(parentDir,ptID));
    %% [20181204] initializing important variables and looping over those files.
    if ~exist(ptWaveFile,'file')
        % initializing data matrices for wave velocities and directions.
        patientWaveBetas = cell(1,length(dirList)-2);
        patientWaveDirections = cell(1,length(dirList)-2); %-2 for ./ and ../
        patientWaveSpeeds = cell(1,length(dirList)-2);
        patientWaveTimes = cell(1,length(dirList)-2);
        patientWavePs = cell(1,length(dirList)-2);
        patientWaveNPDirections = cell(1,length(dirList)-2); %-2 for ./ and ../
        patientWaveNPSpeeds = cell(1,length(dirList)-2);
        
        
        % [20181017] looping over files in each IED directory
        dirTic = tic;
        for fl = 3:length(dirList)
            % loading preprocessed IED data.
            fprintf('\nloading data file %d of %d for patient %s...',fl-2,length(dirList)-2,ptID)
            waveDataTic = tic;
            load(fullfile(dirList(fl).folder,dirList(fl).name));
            A = toc(waveDataTic);
            fprintf('\n     ...took %.2f seconds.',A)
            
            %% [20181119] recovernig code from IEDwaves in order to recover the detection timepoints.
            % see IEDwaves_U.m for more info.
            detectionThreshold = 5;
            nSamps = IEDdata.parameters.downSamplingRate/2; % note this is double the previous window size (this is 1/2 second)
            nBins = ceil(IEDdata.resampledDataLength/nSamps);
            allDetections = sort(cat(1,IEDdata.detections.times));
            [detectionHisto,detectionEdges] = histcounts(allDetections,nBins);
            relevantLefts = detectionEdges([detectionHisto false]>detectionThreshold);
            relevantRights = detectionEdges([false detectionHisto]>detectionThreshold);
            
            % now find the detections within these bin limits.
            retainedDetectionIdcs = (allDetections>=relevantLefts & allDetections<=relevantRights);
            
            
            %% [20181119] if there is wave data, find speed, direction and timing, and save it to new variables for each patient.
            if isfield(IEDdata,'waves')
                % initializing data matrices for wave velocities and directions.
                patientWaveBetas{fl-2} = nan(3,length(IEDdata.waves));
                patientWaveDirections{fl-2} = nan(1,length(IEDdata.waves)); %-2 for ./ and ../
                patientWaveSpeeds{fl-2} = nan(1,length(IEDdata.waves));
                patientWaveTimes{fl-2} = nan(1,length(IEDdata.waves));
                patientWavePs{fl-2} = nan(1,length(IEDdata.waves));
                patientWaveNPDirections{fl-2} = nan(1,length(IEDdata.waves));
                patientWaveNPSpeeds{fl-2} = nan(1,length(IEDdata.waves));
                
                % getting directions of all significant and non-artifactual waves
                for wv = 1:length(IEDdata.waves)
                    if (~isempty(IEDdata.waves(wv).p))  % && lt(IEDdata.waves(wv).p,0.05) % need a way of determinnig which are actual IEDs, and which are artifacts
                        patientWaveBetas{fl-2}(:,wv) = pinv(IEDdata.waves(wv).beta);
                        IEDdata.waves(wv).V = pinv(IEDdata.waves(wv).beta(1:2));
                        patientWaveDirections{fl-2}(wv) = IEDdata.waves(wv).V(1) + IEDdata.waves(wv).V(2)*sqrt(-1); % converting to complex numbers
                        patientWaveSpeeds{fl-2}(wv) = sqrt(sum(IEDdata.waves(wv).V.^2));
                        patientWaveTimes{fl-2}(wv) = median(allDetections(retainedDetectionIdcs(:,wv)));
                        patientWavePs{fl-2}(wv) = IEDdata.waves(wv).p;
                        patientWaveNPDirections{fl-2}(wv) = IEDdata.nonParamWaves(wv).direction;
                        patientWaveNPSpeeds{fl-2}(wv) = IEDdata.nonParamWaves(wv).speed;
                    end
                end
            else
                fprintf('\nno wave data for file %d of %d for patient %s...',fl-2,length(dirList)-2,ptID)
            end
        end
        B = toc(dirTic);
        fprintf('\noverall traveling wave parsing took %.2f minutes.',B/60)
        
        
        %% saving data.
        save(ptWaveFile,'patientWaveBetas','patientWaveDirections','patientWaveSpeeds','patientWaveTimes','patientWavePs','patientWaveNPDirections','patientWaveNPSpeeds','-V7.3')
    else
        % loading already computed data [saves a lot of time]
        load(ptWaveFile)
    end
    
    
    %% getting the directions of seizure expansion from the ictal wavefront
    if ptHasSz
        szFileDir = '/home/user1/data/Seizures/UofU';
        for sz = 1:length(seizureFiles)
            load(fullfile(szFileDir,ptID,seizureFiles{sz}))
            
            % generating smoothed (500 ms) psths
            % clearing variables
            clear MUAdata rasterSpikes tmpPSTHs tPSTH channelPSTHs tmp
            for ch = 1:length(goodChannels)
                MUAdata.times = mua_data.timestamps{goodChannels(ch)};
                rasterSpikes{ch} =  MUAdata.times;
                [tmpPSTHs(goodChannels(ch),:),tPSTH,~] = psth(MUAdata,.5,'n',[preictalTimes(sz) preictalTimes(sz)+seizureDurations(sz)],1);
                tmp = tPSTH-preictalTimes(sz);
                tPSTH = tmp(tmp>0 & tmp<seizureDurations(sz));
                channelPSTHs(goodChannels(ch),:) = tmpPSTHs(goodChannels(ch),tmp>0 & tmp<seizureDurations(sz));
            end
            channelPSTHs = channelPSTHs(goodChannels,:);
            
            % find maxima during wavefront period.
            [maxAmps,maxIdcs] = max(channelPSTHs(:,tPSTH>IWwins(sz,1) & tPSTH<IWwins(sz,2)),[],2);
            
            %             % find local maxima in firing slopes.
            %             eta = 0.001;
            %             [maxDiffAmps,maxIdcs] = max(gradient(channelPSTHs(:,tPSTH>IWwins(sz,1) & tPSTH<IWwins(sz,2))),[],2);
            %
            MAXPLT = false;
            if MAXPLT
                % visualizing maximum detection
                figure
                hold on
                plot(tPSTH(tPSTH>IWwins(sz,1) & tPSTH<IWwins(sz,2)),channelPSTHs(:,tPSTH>IWwins(sz,1) & tPSTH<IWwins(sz,2)))
                scatter(tPSTH(maxIdcs)+IWwins(sz,1),maxAmps,10,rgb('black'),'filled')
                %                 scatter(tPSTH(maxIdcs),channelPSTHs(maxIdcs),20,rgb('gold'),'filled')
                hold off
                ylabel('firing rate (spikes/s)')
                xlabel('time relative to seizure onset (s)')
                title('maxima detection')
                keyboard
            end
            
            
            %% wavefront multilinear regression measurements.
            % spatial linear regression on wavefront maxima.
            RESPONSE = num2cell(tPSTH(maxIdcs)+IWwins(sz,1));
            [IWbeta(:,sz), IWV(:,sz), IWp(sz)] = ...
                SpatialLinearRegression(RESPONSE,P,'switch_plot',0,'Lossfun','L1'); %\
            IWV(:,sz) = pinv(IWbeta(1:2,sz));
            IWspeed(sz) = sqrt(sum(IWV(:,sz).^2)); %(in cm/s)
            IWdirection(sz) = angle(IWV(1,sz) + IWV(2,sz)*sqrt(-1));
            
            
            %% wavefront nonparametric methods.
            [~,channelOrder] = sort(tPSTH(maxIdcs),'ascend');    % channelOrder: channels ordered by their burst minimae
            % looping over channels
            for chz = 1:length(channelOrder)
                try
                    [I,J] = ind2sub(size(map),find(map==channelOrder(chz)));
                    % saving subscripts over Channels
                    Ich(chz) = I;
                    Jch(chz) = J;
                    ArrayOrder(I,J) = chz;
                catch
                    fprintf('couldn"t find channel %d in the map.',chz);
                end
            end
            
            
            %% calculating burst direction.
            binaryImage = ArrayOrder~=0;
            measurementsPlus = regionprops(logical(binaryImage),ArrayOrder,'WeightedCentroid');
            measurementsMinus = regionprops(logical(binaryImage),length(channelOrder)-ArrayOrder,'WeightedCentroid');
            
            Rminus = measurementsMinus.WeightedCentroid;
            Rplus = measurementsPlus.WeightedCentroid;
            
            propDirection = Rplus-Rminus;
            
            [NPIWdirection(sz),NPIWspeed(sz)] = cart2pol(propDirection(1),propDirection(2));
            
            % visualizing seizure
            pltSz = false;
            if pltSz
                figure(pt*1000+sz)
                plotmultipleaxes(1,2,3,0.05,pt*1000+sz)
                hold on
                imagesc(tPSTH,goodChannels,channelPSTHs)
                line([IWwins(sz,1) IWwins(sz,1)],[1 max(goodChannels)],'linestyle','--','color',rgb('silver'))
                line([IWwins(sz,2) IWwins(sz,2)],[1 max(goodChannels)],'linestyle','--','color',rgb('silver'))
                scatter(IWwins(sz,1)+tPSTH(maxIdcs),1:length(maxIdcs),2,rgb('white'),'filled')
                hold off
                xlabel('time relative to seizure onset (s)')
                ylabel('microelectrode channels')
                axis tight square
                title(sprintf('neuronal firing, %s seizure: %d',ptID,sz))
                
                plotmultipleaxes(2,2,3,0.05,pt*1000+sz)
                plotSpikeRaster(rasterSpikes,'PlotType','vertline');
                xlim([preictalTimes(sz) preictalTimes(sz)+seizureDurations(sz)])
                axis off square
                
                plotmultipleaxes(3,2,3,0.05,pt*1000+sz)
                hold on
                patch([tPSTH fliplr(tPSTH)],[mean(channelPSTHs)+std(channelPSTHs)./sqrt(length(goodChannels)) fliplr(mean(channelPSTHs)-std(channelPSTHs)./sqrt(length(goodChannels)))],rgb('black'),'edgecolor','none','facealpha',0.5)
                plot(tPSTH,mean(channelPSTHs),'color',rgb('black'))
                hold off
                xlabel('time relative to seizure onset (s)')
                ylabel('firing rate (spikes/s)')
                axis square
                
                plotmultipleaxes(3,2,2,0.1,pt*1000+sz)
                histogram(tPSTH(maxIdcs),20)
                xlabel('timing of maximum firing rate')
                ylabel('count')
                axis square
                title('histogram of delays across the array')
                
                plotmultipleaxes(4,2,2,0.1,pt*1000+sz)
                scatter3(P(:,1),P(:,2),cell2mat(RESPONSE(:)),'filled');hold on;
                [P1U,P2U] = meshgrid(sort(unique(P(:,1))),sort(unique(P(:,2))));
                f = scatteredInterpolant(P(:,1),P(:,2),P*IWbeta(1:2,sz) + IWbeta(end,sz));
                Z = f(P1U,P2U);
                mesh(P1U,P2U,Z) %interpolated
                xlabel('cm');ylabel('cm');zlabel('Second');
                title(sprintf('speed: %.2f cm/s, p = %.2f',IWspeed(sz),IWp(sz)))
                
                % saving
                halfMaximize(pt*1000+sz,'left')
                print(pt*1000+sz,sprintf('/media/user1/data4TB/Figs/Seizures/%s/%s_ictalWavefrontDetails_seizure%d',ptID,ptID,sz),'-dpdf','-fillpage')
            end
            
        end
    end
    
    
    %% figures within patients [quantify speed and direction w/in pt]
    waveBetas = cell2mat(patientWaveBetas);
    wavePs = cell2mat(patientWavePs);
    waveDirs = cell2mat(patientWaveDirections);
    waveSpds = cell2mat(patientWaveSpeeds);
    isWave = wavePs<0.05; % alpha = 0.05
    waveDirsNP = cell2mat(patientWaveNPDirections);
    waveSpdsNP = cell2mat(patientWaveNPSpeeds);
    
    % removing outliers.
    idx = outliers(waveSpds,[25 75],3);
    isWave(idx) = false;
    isWave(waveSpds>200) = false;
    
    
    nonparamPLT = false;
    %% [20181120] it turns out that eta doesn't have any traveling waves.
    if ~isempty(waveDirs)
        figure(pt)
        % compass plotting wave directions
        subplot(2,2,1)
        cH = compass(waveBetas(1,isWave) + waveBetas(2,isWave)*sqrt(-1));
        set(cH,{'Color'},num2cell(copper(sum(isWave)),2))
        xlabel('IED traveling wave directions')
        title('multilinear regression model')
        
        % histogram of wave speeds
        subplot(2,2,2)
        hold on
        hH = histfit(waveSpds(isWave),10,'kernel');
        tH = text(mean(waveSpds(isWave)),15,sprintf('mean+/-std = %.2f+/-%.2f',mean(waveSpds(isWave)),std(waveSpds(isWave))));
        % deets
        hH(1).FaceColor = rgb('darkslategray');
        hH(2).Color = rgb('darkkhaki');
        tH.FontWeight = 'bold';
        tH.Color = rgb('darkkhaki');
        ylabel('count')
        xlabel('traveling wave speeds (cm/s)')
        axis square
        title('multilinear regression model')
        
        % polar histogram of wave directions.
        subplot(2,2,3)
        polarhistogram(angle(waveDirs(isWave)),18,'EdgeColor',rgb('darkslategray'),'FaceColor',rgb('darkslategray'),'FaceAlpha',0.5)
        title('circular histogram of wave angles (20 degrees/bin)')

        % waveFRONT directions
        subplot(2,2,4)
        try
            iwH = compass(IWV(1,:) + IWV(2,:)*sqrt(-1));
            xlabel('wavefront expansion direction')
            title('ictal wavefront directions')
        catch
            title('no seizures for this patient... yet')
        end
        
        %saving figures for each patient
        halfMaximize(pt,'left')
        saveDir = sprintf('/media/user1/data4TB/Figs/IEDs/%s',ptID);
        figName = sprintf('%s_allIEDwaveSpeedsAndDirections.pdf',ptID);
        print(fullfile(saveDir,figName),'-dpdf','-fillpage')
        
        %% also saving to my DB
        print(fullfile('~/Dropbox',figName),'-dpdf','-fillpage')
        
        
        %% [20181206] for nonparametric wave measures.
        if nonparamPLT
            figure
            [waveCartNPX, waveCartNPY] = pol2cart(waveDirsNP(isWave),waveSpdsNP(isWave));
            
            subplot(2,2,1)
            cHNP = compass(waveCartNPX, waveCartNPY);
            set(cHNP,{'Color'},num2cell(redgrayblue(sum(isWave)),2))
            xlabel('IED traveling wave directions')
            title('nonparametric model')
            
            % plotting wave speeds
            subplot(2,2,2)
            hold on
            hH = histfit(waveSpdsNP(isWave),10,'kernel');
            tH = text(mean(waveSpdsNP(isWave)),max(hH(2).YData)+2,sprintf('mean+/-std = %.2f+/-%.2f',mean(waveSpdsNP(isWave)),std(waveSpdsNP(isWave))));
            % deets
            hH(1).FaceColor = rgb('rosybrown');
            hH(2).Color = rgb('orangered');
            tH.FontWeight = 'bold';
            tH.Color = rgb('orangered');
            ylabel('count')
            xlabel('traveling wave speeds (relative units)')
            axis square
            title('nonparametric model')
            
            subplot(2,2,3)
            polarhistogram(waveDirsNP(isWave),18,'EdgeColor',rgb('darkslategray'),'FaceColor',rgb('darkslategray'),'FaceAlpha',0.5)
            
            % waveFRONT directions
            subplot(2,2,4)
            try
                [iwX,iwY] = pol2cart(NPIWspeed,NPIWdirection);
                iwH = compass(iwX,iwY);
                xlabel('wavefront expansion direction')
                title('ictal wavefront directions')
            catch
                title('no seizures for this patient... yet')
            end
            
            %saving figures for each patient
            halfMaximize(pt,'left')
            saveDir = sprintf('/media/user1/data4TB/Figs/IEDs/%s',ptID);
            figName = sprintf('%s_allIEDwaveSpeedsAndDirections_nonparametricMeasures.pdf',ptID);
            print(fullfile(saveDir,figName),'-dpdf','-fillpage')
        end
    end
    %% [20181206] statistics
    
end


%% figures across patients


%% [20181105] need to figure out how inlcuding the unit data will improve the story...



