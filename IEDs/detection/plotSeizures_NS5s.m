% script to open and plot patient UT array data if needed
clear; close all;

nPts = 6;
for pt = 3 %2:nPts
    % (s)which patients
    clearvars -except parentDir pt nPts
    switch pt
        case 1
            ptID = 'nu';
            % sad magic variables resulting from my lack of prescience go here
            Fs = 3e4;
        case 2
            ptID = 'gamma';
            Fs = 3e4;
        case 3
            ptID = 'Epsilon';
            Fs = 3e4;
        case 4
            ptID = 'eta';
            Fs = 3e4;
        case 5
            ptID = 'kappa';
            Fs = 3e4;
        case 6
            ptID = 'zeta';
            Fs = 3e4;
    end
    
    % where to save figs.
    patientDir = '/home/user1/Desktop/R5/Tyler/Seizure/';
    
    % finding all of the files
    if strcmp(patientDir(end),'/')
        files = subdir([patientDir ptID '/*.ns5']);
    else
        files = subdir([patientDir '/' ptID '/*.ns5']);
    end
    
    tmpSzDir = '~/Dropbox/USzPics/';
    % looping over files
    for fl = 1:length(files)
        fprintf('\nopening file: %s',files(fl).name)
        nSegs = 8;
        NS5 = openNSx(files(fl).name,'read','s:3e4');
        
        % setting up axes
        nChans = size(NS5.Data,1);
        
        dataSize = size(NS5.Data,2)*3e4;
        chunkSize= floor(dataSize./nSegs);
        for chk = 1:nSegs
            fprintf('\nloading data from file %s \n...for chunk %d of %d...\n',files(fl).name,chk,nSegs)
            
            % opening a ~ 8 - 10 GB chunk of data
            NS5 = openNSx(files(fl).name,'read',['t:' num2str(1+chunkSize*(chk-1)) ':' num2str(chunkSize*chk)],'sample');
            
            % plotting the 
            tSec = linspace(0,size(NS5.Data,2)./3e4,size(NS5.Data,2));
            
            % plotting images of all channels.
            figure
            imagesc(tSec,1:nChans,NS5.Data,[-2e3 2e3])
            colorbar
            ylabel('channels')
            xlabel('time (s)')
            
            % getting current file name
            [fPath,fName,ext] = fileparts(files(fl).name);
            
            % saving pictures.
            halfMaximize(gcf,'left')
            print(gcf,sprintf('%s%s_dataFromFile_%s_chunk%d_of%d',tmpSzDir,ptID,fName,chk,nSegs),'-dpng')
            close(gcf)
        end
    end
    
end