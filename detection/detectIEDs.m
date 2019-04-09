function [IEDdata] = detectIEDs(data,Fo)
% DETECTIEDS finds interictal epileptiform discharges (IEDs)
%
%   detectIEDs(data) returns the times of IEDs found in the matrix [data::
%   channels X samples]
%

% Author: [EHS::20181018]

IEDdata.parameters.OGdataLength = length(data);
IEDdata.parameters.OGsamplingRate = Fo;

% detection parameters
Fs = 400;
IEDdata.parameters.downSamplingRate = Fs;
IEDdata.parameters.detectionThreshold = 8;
IEDdata.parameters.artifactThreshold = 20;
IEDdata.parameters.windowSize = 1;      % full window around each IED in seconds.

% filtering and downsampling
[b,a] = butter(4,[20 40]/(Fs/2));
[nChans nSamps] = size(data);
% making a binary matrix for IED detections.
IEDmat = false(nChans,floor(nSamps./(Fo/Fs)));
tSec = linspace(0,nSamps./Fs,floor(nSamps./(Fo/Fs)));
for ch = 1:nChans
    % filtering and thresholding
    updateUser(ch,10,nChans,'downsampling and filtering channel')
    tmpSig = resample(data(ch,:),Fs,Fo)-mean(resample(data(ch,:),Fs,Fo));
    data2040 = filtfilt(b,a,tmpSig);
    IEDdata.resampledDataLength = length(tmpSig);
    IEDdata.measurements(ch).SNR = IEDdata.parameters.detectionThreshold*std(abs(data2040));
    IEDdata.measurements(ch).ArtifactThresh = IEDdata.parameters.artifactThreshold*std(abs(data2040));
    
    % finding peaks
    [IEDdata.foundPeaks(ch).peaks,IEDdata.foundPeaks(ch).locs,IEDdata.foundPeaks(ch).peakWidth,IEDdata.foundPeaks(ch).peakProminence]...
        = findpeaks(smooth(abs(data2040),Fs/50),'MinPeakHeight',IEDdata.measurements(ch).SNR);
    
    % [20180925] retain those detections occurring within 250 ms
    locInds = ~([0; diff(IEDdata.foundPeaks(ch).locs)]<(IEDdata.parameters.windowSize*Fs/4));
    
    % retaining time stamps.
    chosenLocs = IEDdata.foundPeaks(ch).locs(locInds);
    locSize(ch) = length(chosenLocs);
    IEDdata.detections(ch).times = chosenLocs;
    
    % clipping data around these peaks.
    nPeaks(ch) = length(chosenLocs);
%     if gt(nPeaks(ch),max(nPeaks(1:ch-1)))
        for pk = 1:nPeaks(ch)
            updateUser('extracting peak',pk,10,nPeaks(ch))
            if and(gt(chosenLocs(pk)-floor(IEDdata.parameters.windowSize*Fs/2),0),lt(chosenLocs(pk)+ceil(IEDdata.parameters.windowSize*Fs/2),length(tmpSig)))
                % storing all channels of full bandwidth data. wish me luck
                IEDdata.fullBWdata(ch,pk).windowedData = ...
                    data(:,(Fo/Fs)*chosenLocs(pk)-floor(IEDdata.parameters.windowSize*Fo/2):...
                    (Fo/Fs)*chosenLocs(pk)+ceil(IEDdata.parameters.windowSize*Fo/2));
                IEDdata.resampledData(ch,pk).windowedData = ...
                    tmpSig(chosenLocs(pk)-floor(IEDdata.parameters.windowSize*Fs/2):...
                    chosenLocs(pk)+ceil(IEDdata.parameters.windowSize*Fs/2));
            end
            % populating the IED matrix
            IEDmat(ch,chosenLocs(pk)) = true;
            
            % saving time vector.
            IEDdata.tSamps = ...
                chosenLocs(pk)-floor(IEDdata.parameters.windowSize*Fo/2):...
                chosenLocs(pk)+ceil(IEDdata.parameters.windowSize*Fo/2);
            % saving time vector.
            IEDdata.tReSamps = ...
                chosenLocs(pk)-floor(IEDdata.parameters.windowSize*Fs/2):...
                chosenLocs(pk)+ceil(IEDdata.parameters.windowSize*Fs/2);

        end
%     else
%         fprintf('\nfewer peaks detected on channel %d, so not saving these data',ch)
%     end
    clear chosenLocs
end

 

%% [20181018] visualization stuff below here. Would be nice to get this going for the full detection run. 

% % generalTiming
% if exist('pk','var')
%     tSec = linspace(0,IEDdata.parameters.windowSize,length(IEDdata.tSamps));
% end
% 
% % picking the median channel.
% [~,medChans] = sort(locSize);
% 
% keyboard
% 
% %% visualizing IEDs
% % now useless! 
% figure
% plotmultipleaxes(1,1,2,0.07,gcf)
% plotSpikeRaster(IEDmat,'PlotType','vertline');
% axis square tight
% xlabel('time (samples)')
% ylabel('channels')
% plotmultipleaxes(2,1,2,0.07,gcf)
% for pk2 = 1:length(IEDdata.fullBWdata)
%     plot(tSec,mean(IEDdata.fullBWdata(pk2).peakData))
% end
% title(['detection S.D. = ' num2str(IEDdata.parameters.detectionThreshold)])
% axis square tight
% xlabel('time (samples)')
% ylabel('LFP (uV)')
% 
% halfMaximize(gcf,'right')




