function IEDtimePlot(ptID)
% IEDtimePlot will dig into a directory of data from the IED analysis
% pipeline and plot various features of the IED data over the full duration
% of the data. 

% author:: EHS20190220

parentDir = '~/data/IEDs';
set(0,'defaultfigurerenderer','painters')


%% [20181017] to start: dig through the IED detections. (separate by patient?)
nPts = 6;
for pt = 2
    % (s)which patients
    clearvars -except parentDir pt nPts 
    switch pt
        case 1
            ptID = 'nu';
            % sad magic variables resulting from my lack of prescience go here
            Fs = 3e4;
            negativeGoingIED = true;
        case 2
            ptID = 'gamma';
            Fs = 3e4;
            negativeGoingIED = true;
        case 3
            ptID = 'epsilon';
            Fs = 3e4;
            negativeGoingIED = true;
        case 4
            ptID = 'eta';
            Fs = 3e4;
            negativeGoingIED = true;
        case 5
            ptID = 'kappa';
            Fs = 3e4;
            negativeGoingIED = false;
        case 6
            ptID = 'zeta';
            Fs = 3e4;
            negativeGoingIED = true;
    end


        % [20181017] looping over files in each IED directory
    dirList = dir(fullfile(parentDir,ptID));
    for fl = 3:length(dirList)
        % loading preprocessed IED data.
        load(fullfile(dirList(fl).folder,dirList(fl).name));

        % getting timing of the start of the file. 
        fileTimeStr = dirList(fl).name(12:end-4);
        
        keyboard
    end
end


        
