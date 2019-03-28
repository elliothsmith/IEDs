% digging through IED detections to start quantifying traveling waves.
%% This script is fast if not calculating waves. Slower if recaluclating vars.
% how slow: 20 mins to run the first three pateints, calculating
% everyhting.

%% [20181206] make sure you've used the distribution version of Jyunyou's
%% multilinear regression code to calculate traveling wave betas and velocity.

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
        [Ch,IA,IB] = intersect(map,goodChannels);
        [p1,p2] = ind2sub([10 10],IA);
        P = [p1,p2]*0.04; % cm
    else
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
        
        
        %% determining the relative timing of every event. 
        % first pulling out the start time from the file name.
        startDate = dirList(fl).name(13:20);
        startTime = dirList(fl).name(22:27);
        detectionTimes = allDetections.IEDdata.parameters.downSamplingRate;
        
        
        
        
        
        
        %% [20181119] if there is wave data, find speed, direction and timing, and save it to new variables for each patient.
        if isfield(IEDdata,'waves')
            % getting directions of all significant and non-artifactual waves
            for wv = 1:length(IEDdata.waves)
                if (~isempty(IEDdata.waves(wv).p))  % && lt(IEDdata.waves(wv).p,0.05) % need a way of determinnig which are actual IEDs, and which are artifact
                    
                    % dataMap(p1(ic),p2(ic)) = IEDdata.fullBWdata(wv).windowedData(ic,sampRange(ff))
                    
                end
            end
        end
    end
end
