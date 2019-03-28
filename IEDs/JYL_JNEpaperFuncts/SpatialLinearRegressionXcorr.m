function [ Beta, Sigma, pValue ] = SpatialLinearRegressionXcorr( LFP, P, varargin )
% [ Rec ] = SpatialLinearRegressionXcorr( LFP, P, varargin )
% 
% Using cross-correlation method to calculate traveling wave velocity
%
% Input: LFP, an n_channel by n_time numeric matrix 
%        P, a n_channel by 2 matrix represents each channel's position
%        [varargin] - 'weighted': true / false, true: it will be weighted 
%                     'Lossfun' : L1 / L2
%                     'n_shuffle': 0 or empty - use F-test, otherwise, use
%                                  shuffle test
%                     'switch_plot': true / false
%                     'maxlag': default: 0.05 sec
%                     'Fs': sampling rate
%
% Output: Beta, a n_parameter by 1 vector
%         Sigma, estimator covariance matrix
%         pValue, f-test result

p = inputParser;
p.CaseSensitive = false;
addParameter(p,'weighted',false);
addParameter(p,'switch_plot',0); 
addParameter(p,'n_shuffle',100); 
addParameter(p,'Lossfun','L2'); 
addParameter(p,'maxlag',0.05); % unit: second
addParameter(p,'Fs',1000); 
parse(p,varargin{:});

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

Sel = ~isnan(dX1);
dX = [dX1(Sel) dX2(Sel)]; % unit: cm

% Regression
dT = Lag(Sel)/p.Results.Fs; % Unit: second
W = Cmax(Sel);
W(W<0) = 0;

if p.Results.switch_plot
    f = figure;A=axes;
    scatter3(dX(:,1),dX(:,2),dT,20);
end

if p.Results.weighted 
    dX = bsxfun(@times,dX,sqrt(W));
    dT = bsxfun(@times,dT,sqrt(W));
end

if strcmpi(p.Results.Lossfun,'L2')
    [Beta,~,~,Sigma] = lscov(dX,dT);
    % F-test
    [Xm,Xn] = size(dX);
    F_CDF = @(f,df1,df2) betainc(f/(f+df2/df1),df1/2,df2/2);    
    Fstat = ((var(dT,0)-var(dT-dX*Beta,0))/Xn) ...
           /(var(dT-dX*Beta,0)/(Xm-Xn)); % unweighted, no more -1 because we do not estimate the intercept
    try
        pValue = 1- F_CDF(Fstat,Xn,Xm-Xn); % F-test
    catch
        pValue = 1;
    end
elseif strcmpi(p.Results.Lossfun,'L1')
    [Beta, ~, ~, pValue, Sigma] = L1LinearRegression(dX,dT);
end

if p.Results.switch_plot
    mesh(A,dX(:,1),dX(:,2),dX*Beta);
end

end


%% Local function - L1 regression
function [B, E, E_s, L1_p_shuffle, Sigma] = L1LinearRegression(X,Y,Weight) 
    if nargin < 3
        Weight = 1;
    end
    % Determine size of predictor data 
    [n,m] = size(X); 
    % Here you should not add the constant term
    X = [X];
    [B,E,Sigma] = L1Optimize(X,Y,Weight);
    p.Results.n_shuffle = 200;
    if p.Results.n_shuffle
        % Shuffle the channel's spatial information
        E_s = nan(p.Results.n_shuffle,1);
            for i_shuffle = 1:p.Results.n_shuffle                 
                [~,E_s(i_shuffle)] = L1Optimize(Shuffle_Position(X),Y,Weight);
            end
            L1_p_shuffle = sum(E_s<=E)/p.Results.n_shuffle;
    else
        L1_p_shuffle = nan;

    end        
    function [Bout,Eout, Sigma] = L1Optimize(Xin,Yin,Win)
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
            Bout = (repmat(w,[1 m]) .* Xin) \ (w .* Yin); % not m+1 because there is no constant term
            iter = iter + 1;
            if iter > 30;break;end
        end 
        % Report mean absolute deviation
        Eout = sum(w.* abs(Xin * BOld - Yin));
        % Calculate covariance
        [~,~,~,Sigma] = lscov(Xin,Yin,w);
    end
end

function Xout = Shuffle_Position(Xin)
    [X_u,~,IB_X] = unique(Xin,'rows');
    Per = randperm(size(X_u,1));
    X_u_per = X_u(Per,:);
    Xout = X_u_per(IB_X,:);                          
end