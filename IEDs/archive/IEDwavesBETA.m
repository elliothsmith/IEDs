% digging through IED detections to start quantifying traveling waves.

parentDir = '~/data/IEDs';
    
%% [20181017] to start: dig through the IED detections. (separate by patient?)
nPts = 2;
for pt = 1:nPts
    % (s)which patients
    clearvars -except parentDir pt nPts
    switch pt
        case 1
            ptID = 'nu';
            % sad magic variables resulting from my lack of prescience go here
            Fs = 400;
        case 2
            ptID = 'gamma';
            Fs = 400;
    end

    % looping over files in each IED directory
    dirList = dir(fullfile(parentDir,ptID));
    T = cell(96,1);
    for fl = 3:length(dirList)
        % loading data
        load(fullfile(dirList(fl).folder,dirList(fl).name));
        
        % find detections that occur over a certain number of channels.
        % recreate the IED mat to simply sum over rows
        % [20181017] not the most efficient...
        IEDmat = false(length,IEDdata.resampledDataLength);
        for ch = 1:length(IEDdata.IEDmat)
            if ~isempty(IEDdata.IEDmat(ch).times)
                for pk = 1:length(IEDdata.IEDmat(ch).times)
                    % populating the IED matrix
                    IEDmat(ch,IEDdata.IEDmat(ch).times(pk)) = true;
                end
            end
        end
        
        % visualize detections again.
        visualizeDetections = false;
        if visualizeDetections
            figure
            % raster plot of IED detections
            plotmultipleaxes(1,1,2,0.07,gcf)
            plotSpikeRaster(IEDmat,'PlotType','vertline');
            axis square tight
            xlabel('time (samples)')
            ylabel('channels')
            % time histogram of detections at millisecond resolution (TOO TIGHT!)
            plotmultipleaxes(2,1,2,0.07,gcf)
            plot(sum(IEDmat))
            axis square tight
            xlabel('time (samples)')
            ylabel('LFP (uV)')
            halfMaximize(gcf,'right')
        end
        
        keyboard
        % detecting timepoints with simultaneously-detected discharges
        chanThreshold = 2; % starting very low for the first one.
        wavePoints = find(sum(IEDmat)>chanThreshold);
        
        nWavePoints = length(wavePoints);
        for  wp = 1:nWavePoints
            
            % detect maximal descent
            maxDescents = min(diff(IEDdata.IEDmat(wavePoints(wp)-winSize*Fs/2)));
            
            % % build timing cell array here
            % T{ch} =

            % calculate traveling waves using Negative peak & maximal descent estimator
            [beta, V, p] = SpatialLinearRegression(T,P,'switch_plot',1,'Lossfun','L1') % Least absolute deviation estimator
            
        end
        
        % code form JYL's sample script.
        %         Beta = [0.5;1.5]; % this is the ground truth for beta
        %         T = P * Beta; % This is the timing of negative peak/maximal descent
        %         T = T + randn(size(T)); % Get some noise
        %         T = mat2cell(T,ones(size(T,1),1),size(T,2)); % Put those spikes in cells
        %[beta, V, p] = SpatialLinearRegression(T,P,'switch_plot',1) % Least square estimator
        
        
        % [20181017] NOTE:: There are several items I may have to add to
        % may detection scripts if I want to realign the data or look at
        % multiunits. The other option for doing this, is jsut going back
        % to the original data, which may even be more efficient. Will have
        % to thinkabout this a bit.
        
    end
    
    
    
    
end
