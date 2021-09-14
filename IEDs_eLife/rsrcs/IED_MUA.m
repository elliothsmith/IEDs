function [mua,timestamps,thresh,wfs] = IED_MUA(data,wfFlag)

% [mua,timestamps,thresh,wfs] = IED_MUA(data) implements a very simple multiunit action potential detection
%
% [mua,timestamps,thresh,wfs] = IED_MUA(data,wfFlag) saves detected waveforms if wfFlag is true (default: false)
%	the output 'wfs' variable has one cell per channel and the data in each cell is a matrix of the form [samples(n=60)  X spike detections]

if isequal(nargin,1)
	wfFlag = false; 
	wfs = [];
end

Fs = 30000;
% number of channels
numChans=size(data,1);

% Band pass between 500 and 3000 Hz
band = [300 3000];
[b,a] = fir1(96,[band(1)/(Fs/2) band(2)/(Fs/2)]);
for ch3=1:numChans
	mua(ch3,:) = filtfilt(b,a,double(data(ch3,:)));
 	% display(['filtering channel ' num2str(ch3)])

	% removing the mean from each channel's signal.
	mua(ch3,:) = mua(ch3,:)-mean(mua(ch3,:));

	% finding a threshold based on the period defined in ThreshPeriod
	% thresh = -4*median(abs(mua(ch3,:)./0.6745));
	thresh = -3.5*median(abs(mua(ch3,:)./0.6745));
	% thresh = -3.5*rms(mua(ch3,:));

	% finding the minima of the filtered signal.
	mua_peaks = find_inflections(mua(ch3,:),'minima');
	fprintf('filtered and thresholded channel %d. Threshold value (uV): %.2f', num2str(ch3),thresh,ch3)

	% saving the times of the minima that are greater than the threshold.
	spiketimes = mua_peaks(mua(ch3,mua_peaks)<thresh);
	timestamps{ch3} = spiketimes./Fs;

	if wfFlag
			wfs = cell(1,length(spiketimes));
		for sp = length(spiketimes):-1:1
			try
				tmp(:,sp) = mua(ch3,spiketimes(sp)-20:spiketimes(sp)+40);
				wfs{ch3} = tmp;
				clear tmp
			end
		end
	end
end

end
