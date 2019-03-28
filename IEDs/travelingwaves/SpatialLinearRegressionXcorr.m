function [ Beta, V, pValue ] = SpatialLinearRegressionXcorr( LFP, P, varargin )
%  [ Beta, V, pValue ] = SpatialLinearRegressionXcorr( LFP, P, varargin )
% 
% Estimating traveling wave velocity by cross-correlation method.
%
%
% Input: 1) LFP, an n_channel by n_time numeric matrix 
% 
%        2) P: The position information. Format: an n_channel by 2 numetric
%           matrix.  For example, a 2 by 2 grid will be presented as
%
%           P = [1 1;
%                1 2;
%                2 1;
%                2 2];
%        
%        3) [Optional inputs] 
%           SpatialLinearRegressionXcorr(Input,P,'Parameter_Name',parameter);
%
%           switch_plot: (1/0) - plot the result or not
%           n_shuffle: n of shuffle test 
%           Fs: sampling frequency for time series data - default:1000 Hz
%           Lossfun: Define the loss function in regression
%                    Options: 
%                    'L2' - least squares
%                    'L1' - least absolute deviation
%                    Default: 'L2'
%           maxlag: the search range of cross-correlation, default: 0.05sec
%
% Output: Beta: The estimated regression coefficients 
%         V: The estimated traveling wave velocity
%         pValue: p Value of shuffle test
%
% Author: Jyun-you Liou
% Final update: 2016/11/30 version 1

p = inputParser;
p.CaseSensitive = false;
addParameter(p,'switch_plot',0); 
addParameter(p,'n_shuffle',200); 
addParameter(p,'Lossfun','L2'); 
addParameter(p,'maxlag',0.05); % unit: second
addParameter(p,'Fs',1000); 
parse(p,varargin{:});

if rank([ones(size(P,1),1), P]) <= size(P,2)
    error('The spatial information contains less dimensions than it shows.')
end

%% Get physical location information
n_channel = size(LFP,1);

%% Episode by episode get cross-correlation matrix
Cmax = nan(n_channel);
Lag = nan(n_channel);
dX1 = nan(n_channel);
dX2 = nan(n_channel);
% Generate pairs
for m1 = 1:n_channel
    for m2 = 1:n_channel
        if m1 >= m2;continue;end
        [C,lag] = xcorr(LFP(m1,:),LFP(m2,:),round(p.Results.maxlag * p.Results.Fs),'coeff');
        [Cmax(m1,m2),Cmaxp] = max(C);
        Lag(m1,m2) = -lag(Cmaxp);
        dX1(m1,m2) = P(m2,1) - P(m1,1);
        dX2(m1,m2) = P(m2,2) - P(m1,2);
    end
end

% The estimation
Sel = ~isnan(dX1);
dX = [dX1(Sel) dX2(Sel)]; % unit: cm

% Regression
dT = Lag(Sel)/p.Results.Fs; % Unit: second

if p.Results.switch_plot
    f = figure;A=axes;
    scatter3(dX(:,1),dX(:,2),dT,20);hold on;
end

%% Since this data is correlated, you need to do permutation test ...
if strcmpi(p.Results.Lossfun,'L2')
    [Beta,~,E0] = lscov(dX,dT);
elseif strcmpi(p.Results.Lossfun,'L1')
    [Beta,E0] = L1Optimize(dX,dT);
end

if p.Results.switch_plot
    [X1,X2] = meshgrid(dX(:,1),dX(:,2));
    f = @(X1,X2) X1*Beta(1) + X2*Beta(2);
    mesh(A,dX(:,1),dX(:,2),f(X1,X2));zlabel('Time')
end


%% Bootstrap for permutation test
if p.Results.n_shuffle && nargout == 3
    for i = p.Results.n_shuffle:-1:1  
        P_shuffled = Shuffle_Position(P);
        % Generate pairs
        for m1 = 1:n_channel
            for m2 = 1:n_channel
                if m1 >= m2;continue;end
                dX1(m1,m2) = P_shuffled(m2,1) - P_shuffled(m1,1);
                dX2(m1,m2) = P_shuffled(m2,2) - P_shuffled(m1,2);
            end
        end

        % The estimation
        dX = [dX1(Sel) dX2(Sel)]; % unit: cm

        % Regression
        dT = Lag(Sel)/p.Results.Fs; % Unit: second

        if strcmpi(p.Results.Lossfun,'L2')
            [~,~,E(i)] = lscov(dX,dT);
        elseif strcmpi(p.Results.Lossfun,'L1')
            [~,E(i)] = L1Optimize(dX,dT);
        end
    end
    pValue = mean(E0>E);
end
V = pinv(Beta);

end


%% Local function - L1 regression
function [Bout,Eout] = L1Optimize(Xin,Yin,Win)
    if nargin < 3;Win = 1;end
    % Use least-squares fit as initial guess
    Bout = Xin\Yin;             
    % Least squares regression 
    BOld = Bout; 
    BOld(end) = BOld(end) + 1e-5; 
    iter = 1;
    % Force divergence % Repeat until convergence 
    while (max(abs(Bout - BOld)) > 1e-6) 
        % Move old coefficients 
        BOld = Bout; 
        % Calculate new observation weights (based on residuals from old coefficients) 
        w = sqrt(Win ./ max(abs(Xin * BOld - Yin),1e-6));
        % Floor to avoid division by zero 
        % Calculate new coefficients 
        Bout = bsxfun(@times,w,Xin) \ (w .* Yin); % not m+1 because there is no constant term
        iter = iter + 1;
        if iter > 30;break;end
    end 
    % Report mean absolute deviation
    Eout = sum(w.* abs(Xin * BOld - Yin));
end


function Xout = Shuffle_Position(Xin)
    [X_u,~,IB_X] = unique(Xin,'rows');
    Per = randperm(size(X_u,1));
    X_u_per = X_u(Per,:);
    Xout = X_u_per(IB_X,:);                          
end