% Demonstration scripts
% First, let's set the spatial information of the recording array (3 by 3)
P = [1 1;
     1 2
     1 3;
     2 1;
     2 2;
     2 3;
     3 1;
     3 2;
     3 3];

%% Multiunit estimator
Beta = [0.5;1.5]; % this is the ground truth for beta
T = P * Beta;
T = repmat(T,[1,10]); % Each channel can have several events
T = T + randn(size(T)); % Get some noise
T = mat2cell(T,ones(size(T,1),1),size(T,2)); % Put those spikes in cells
T{7} = []; % For example, channel 7 fails to detect any spike
[beta, V, p] = SpatialLinearRegression(T,P,'switch_plot',1) % Least square estimator
[beta, V, p] = SpatialLinearRegression(T,P,'switch_plot',1,'Lossfun','L1') % Least absolute deviation estimator

%% Negative peak & maximal descent estimator
Beta = [0.5;1.5]; % this is the ground truth for beta
T = P * Beta; % This is the timing of negative peak/maximal descent
T = T + randn(size(T)); % Get some noise
T = mat2cell(T,ones(size(T,1),1),size(T,2)); % Put those spikes in cells
[beta, V, p] = SpatialLinearRegression(T,P,'switch_plot',1) % Least square estimator
[beta, V, p] = SpatialLinearRegression(T,P,'switch_plot',1,'Lossfun','L1') % Least absolute deviation estimator


%% High gamma envelope estimator
Beta = [0.5;1.5]; % this is the ground truth for beta
T = P * Beta; % This is the timing of maximal gamma power
% Suppose gamma power is Gaussian distributed around the gamma power peak
T_grid = mean(T) + (-50:50)/10;% Let's look at 101 data points around the ictal discharge, & suppose sampling rate = 10
W = normpdf(bsxfun(@minus,T_grid,T),0,1);
W = W + rand(size(W))*max(W(:))/5; % Give it some background noise
imagesc(W); % This is gamma power
[beta, V, p] = SpatialLinearRegression(W,P,'switch_plot',1,'Fs',10) % Least square estimator
[beta, V, p] = SpatialLinearRegression(W,P,'switch_plot',1,'Lossfun','L1','Fs',10) % Least absolute deviation estimator


%% Cross-correlation estimator
Beta = [0.01;-0.005]; % this is the ground truth for beta
T = P * Beta; % This is the timing of maximal gamma power
% Let make a sample LFP
T_grid = mean(T) + (-100:100)/1000;% Let's look at 201 data points around the ictal discharge, & suppose sampling rate = 1000
W = -normpdf(bsxfun(@minus,T_grid,T),0,0.01);
W = W + randn(size(W))*max(abs(W(:)))/10; % Give it some background noise
plot(W'); % This is LFP (yeah .. very ugly)
[beta, V, p] = SpatialLinearRegressionXcorr(W,P,'switch_plot',1,'Fs',1000) % Least square estimator
[beta, V, p] = SpatialLinearRegressionXcorr(W,P,'switch_plot',1,'Lossfun','L1','Fs',1000) % Least absolute deviation estimator

