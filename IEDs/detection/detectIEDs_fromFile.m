function [IEDdata] = detectIEDs_fromFile(ptID,fileName)
% This function is a wrapper for DETECTIEDS.
% DETECTIEDS operates on a data matrix.
% this function grabs the data from a specific file.
% use the function DETECTIEDS_OEVERDIRS to run this as a worm over a
% specific directory.


% TODO:: remove direct link references in possible.

%% loading data
fDeets = dir(fileName);
if fDeets.bytes<2.5e10 % splitting data into smaller segments if it's too big to read.
    fprintf('\nloading full-bandwidth data...\n')
    NS5 = openNSx(fileName);
    if iscell(NS5.Data)
        if size(NS5.Data{2},1)<50
            fprintf('\nnumber of channels is less than 50. This may not be UMA data...')
        else
            if size(NS5.Data{2},1)>96
                [IEDdata] = detectIEDs(double(NS5.Data{2}(1:96,:)),3e4);
            else
                [IEDdata] = detectIEDs(double(NS5.Data{2}),3e4);
            end
        end
        
    else
        if size(NS5.Data,1)<50
            fprintf('\nnumber of channels is less than 50. This may not be UMA data...')
        else
            if size(NS5.Data,1)>96
                [IEDdata] = detectIEDs(double(NS5.Data(1:96,:)),3e4);
            else
                [IEDdata] = detectIEDs(double(NS5.Data),3e4);
            end
        end
    end
    
    % save figure
    [~,fName,~] = fileparts(fileName);
    %     saveas(gcf,['/media/user1/data4TB/Figs/IEDs/IEDdetectionsFromFile_' fName '.pdf'])
    %     close(gcf)
    
    % saving IED data.
    if exist('IEDdata','var')
        if ~exist(['/media/user1/data4TB/data/IEDs/' ptID '/'],'dir')
            mkdir(['/media/user1/data4TB/data/IEDs/' ptID '/']);
        end
        save(['/media/user1/data4TB/data/IEDs/' ptID '/IEDsfromfile_' fName '.mat'],'IEDdata','-v7.3')
    end
else % split the file into N segments.
    nSegs = 16;
    NS5 = openNSx(fileName,'read','s:3e4');
    if iscell(NS5.Data)
        if size(NS5.Data{2},1)<50
            fprintf('\nnumber of channels is less than 50. This may not be UMA data...')
        else
            nChans = size(NS5.Data{2},1);
            dataSize = size(NS5.Data{2},2)*3e4;
            chunkSize= floor(dataSize./nSegs);
            for chk = 1:nSegs
                fprintf('\nloading data from file %s for chunk %d of %d...\n',fileName,chk,nSegs)
                
                % opening a ~ 8 - 10 GB chunk of data
                NS5 = openNSx(fileName,'read',['t:' num2str(1+chunkSize*(chk-1)) ':' num2str(chunkSize*chk)],'sample');
                
                % detecting IEDs for this chunk
                [IEDdata] = detectIEDs(double(NS5.Data{2}(1:96,:)),3e4);
                
                % save figure
                [~,fName,~] = fileparts(fileName);
                %                 saveas(gcf,['/media/user1/data4TB/Figs/IEDs/IEDdetectionsFromFile_' fName '_chunk_' num2str(chk) '.pdf'])
                %                 close(gcf)
                
                % saving IED data for this chunk
                fprintf('\nsaving data for chunk %d of %d...\n',chk,nSegs)
                if exist('IEDdata','var')
                    save(['/media/user1/data4TB/data/IEDs/' ptID '/IEDsfromfile' fName '__chunk' num2str(chk) '.mat'],'IEDdata','-v7.3')
                end
            end
        end
    else % if the file is huge and the data is not a cell (i.e. has no pauses)
        if size(NS5.Data,1)<50
            fprintf('\nnumber of channels is less than 50. This may not be UMA data...')
        else
            nChans = size(NS5.Data,1);
            dataSize = size(NS5.Data,2)*3e4;
            chunkSize= floor(dataSize./nSegs);
            for chk = 1:nSegs
                fprintf('\nloading data from file %s for chunk %d of %d...\n',fileName,chk,nSegs)
                
                % opening a ~ 8 - 10 GB chunk of data
                NS5 = openNSx(fileName,'read',['t:' num2str(1+chunkSize*(chk-1)) ':' num2str(chunkSize*chk)],'sample');
                
                % detecting IEDs for this chunk
                [IEDdata] = detectIEDs(double(NS5.Data(1:nChans,:)),3e4);
                
                % save figure
                [~,fName,~] = fileparts(fileName);
                %                 saveas(gcf,['/media/user1/data4TB/Figs/IEDs/IEDdetectionsFromFile_' fName '_chunk_' num2str(chk) '.pdf'])
                %                 close(gcf)
                
                % saving IED data for this chunk
                fprintf('\nsaving data for chunk %d of %d...\n',chk,nSegs)
                if exist('IEDdata','var')
                    save(['/media/user1/data4TB/data/IEDs/' ptID '/IEDsfromfile' fName '__chunk' num2str(chk) '.mat'],'IEDdata','-v7.3')
                end
            end
        end
    end
end
