% digging through IED detections to start quantifying traveling waves.
%% This script is fast if not calculating waves. Slower if recaluclating vars.
% how slow: 20 mins to run the first three pateints, calculating
% everyhting.

%% [20181206] make sure you've used the distribution version of Jyunyou's
%% multilinear regression code to calculate traveling wave betas and velocity.

set(0,'defaultfigurerenderer','opengl')

% defining the directory where patientID subdirectories are located.
parentDir = '/media/user1/data4TB/data/IEDs';


%% [20181017] to start: dig through the IED detections. (separate by patient?)
nPts = 5;
IEDwaveCount = zeros(1,nPts); % number of detected IEDs overall.
for pt = 2  % 1:nPts
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
    % where to save figs.
    tmpFigDir = '/media/user1/data4TB/data/IEDs/tmpIEDFigs';
    
    % converting the 2D array map into a 96 X 2 matrix with cm units, P.
    map = electrodepinout; % assuming that all of the utah utah arrays had the same pinout, which is the same as those from the early columbia cases.
    
    % first way of calculating P
    if exist('goodChannels','var')
		nChans = length(goodChannels)
        [Ch,IA,IB] = intersect(map,goodChannels);
        [p1,p2] = ind2sub([10 10],IA);
        P = [p1,p2]*0.04; % cm
    else
		nChans = 96;
        [Ch,IA,IB] = intersect(map,1:96);
        [p1,p2] = ind2sub([10 10],IA);
        P = [p1,p2]*0.04; % cm
    end
    
    
    % accessing saved, preprocessed IED data.
    dirList = dir(fullfile(parentDir,ptID));
    %% [20181204] initializing important variables and looping over those files.
    % [20181017] looping over files in each IED directory
    for fl = 3:length(dirList)
        % loading preprocessed IED data.
        fprintf('\nloading data file %d of %d for patient %s...',fl-2,length(dirList)-2,ptID)
        waveDataTic = tic;
        load(fullfile(dirList(fl).folder,dirList(fl).name));
        A = toc(waveDataTic);
        fprintf('\n     ...took %.2f seconds.',A)
        
        
        %% [20181119] recovernig code from IEDwaves in order to recover the detection timepoints.
        % see IEDwaves_U.m for more info.
        detectionThreshold = 20;
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
            % getting directions of all significant and non-artifactual waves
            for wv = 1:length(IEDdata.waves)
                if (~isempty(IEDdata.waves(wv).p)) && lt(IEDdata.waves(wv).p,0.01) % need a way of determinnig which are actual IEDs, and which are artifactsa

					% range of samples to plot
                    sampRange = 1.5e4:1.65e4;

					% First plot examples of the IEDs as a separate figure
					figure(1)
					subplot(1,2,1)
                    plot((1:3e4+1)./3e4,IEDdata.fullBWdata(wv).windowedData(1:nChans,:),'linewidth',0.5) 
					ylabel('LFP (uV)')
					xlabel('time relative to IED detection (s)')
					axis square tight

					subplot(1,2,2)
					plot(sampRange./3e4,IEDdata.fullBWdata(wv).windowedData(1:nChans,sampRange),'linewidth',0.5)
					xlabel('time relative to IED detection (s)')
					axis square tight

					% save
					saveas(1,['~/Dropbox/IEDmovies2/' ptID '_IEDwave' dirList(fl).name 'IEDn' num2str(wv) '.pdf']);
					close(1)

                    %% make a video of the windowed data
                    movie_file_name = ['~/Dropbox/IEDmovies2/' ptID '_IEDwave' dirList(fl).name 'IEDn' num2str(wv) '.avi'];
                    movie_object = VideoWriter(movie_file_name,'Motion JPEG AVI');
                    movie_object.FrameRate = 300;
                    movie_object.Quality = 50;
                    open(movie_object);
                    
                    % setting up figure
                    movie_handle = figure(999);
                    figure(999)
                    maximize(999)
                    
                    % how long is the movie, in samples..
                    cLims = [-500 500];
                    for ff = 1:length(sampRange)
						
						subplot(1,2,1)
						hold on
						plot(mean(IEDdata.fullBWdata(wv).windowedData(:,sampRange)))
						line([sampRange(ff) sampRange(ff)],[min(mean(IEDdata.fullBWdata(wv).windowedData)) max(mean(IEDdata.fullBWdata(wv).windowedData))])
						ylim([min(mean(IEDdata.fullBWdata(wv).windowedData)) max(mean(IEDdata.fullBWdata(wv).windowedData))])
						hold off

						for ic = 1:nChans
                            dataMap(p1(ic),p2(ic)) = IEDdata.fullBWdata(wv).windowedData(ic,sampRange(ff));
                        end
						
						subplot(1,2,2)	
                        surf(interp2(dataMap,5))
						ylim([min(mean(IEDdata.fullBWdata(wv).windowedData)) max(mean(IEDdata.fullBWdata(wv).windowedData))])
                        axis off square

                        % writing video
                        writeVideo(movie_object,getframe(movie_handle));
                        
                    end
                end
            end
        end
    end
end
