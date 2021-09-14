function [Beta, V, pValue] = SpatialLinearRegression(Input,P,varargin)
% [Beta, V, pValue] = SpatialLinearRegression(T,P,varargin)
%
% Estimating traveling wave velocity by multivariate linear regression.
%
% This function is used for the following methods
% 
% - Multiunit estimator
% - Negative peak estimator
% - Maximal descent estimator
% - High gamma envelope estimator
%
% Please see the following description to cater the data format for each
% estimator
%
% Input: 1) Input: 
%           ** Important **
%           Input type can be cell or numeric matrix, depends on which 
%           estimator you are going to use.
%
%           * [Data type: cell] - used for 
%                               * multiunit estimator
%                               * negative peak estimator
%                               * maximal descent estimator
%
%           In this case, Input is 'T' (information of event timing)
%           Input needs to be n_channel by 1 cell, which each cell 
%           corresponds to the information from one recording channel. 
%           Within a cell, it contains a numeric vector which encodes 
%           timing of events (multiunit spike, negative peak, or maximal 
%           descent) happen.
%
%           For example: channel 1 detects 3 multiunit spikes at time 5, 
%                        7, 10 during an ictal discharge.
%                        Input{1} = [5, 7, 10];  
%
%           For example, channel 1 detects its negative peak at time 12
%                        Input{1} = 12;
%
%           * [Data type: numeric matrix] - used for 
%                                         * high gamma envelope estimator
%
%           In this case, Input encodes the instantaneous power. 
%           Input needs to be a n_channel by n_time matrix.  Every entry 
%           of Input needs to be strictly positive.  In other words, Input
%           here is treated as 'W' (weight) for linear regression.
%           
%           For example, channel 1 has instantaneous power ramped up and
%           down within 10 time steps of an ictal discharge
%           
%           Input(1,:) = [1 2 3 4 5 5 4 3 2 1]
%
%           The sampling rate (n_time, the 2nd dimension of Input matrix)
%           is set to be 1000 Hz by default.  You can change it by giving 
%           optional input. For example
%
%           SpatialLinearRegression(Input,P,'Fs',500);
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
%           SpatialLinearRegression(Input,P,'Parameter_Name',parameter);
%
%           switch_plot: (1/0) - plot the result or not
%           n_shuffle: n of shuffle test if shuffle test is called 
%           Fs: sampling frequency for time series data - default:1000 Hz
%           Lossfun: Define the loss function in regression
%                    Options: 
%                    'L2' - least squares
%                    'L1' - least absolute deviation
%                    Default: 'L2'
%
% Output: Beta: The estimated regression coefficients 
%         V: The estimated traveling wave velocity
%         pValue: p Value of F-test or shuffle test
%
% Author: Jyun-you Liou
% Final update: 2016/11/30, version 1

p = inputParser;
addParameter(p,'alpha',0.05); % Significance level
addParameter(p,'switch_plot',0); 
   % In point process, it will plot data average at its physical position.
   % In continuous process, it will return shuffle test results.
addParameter(p,'n_shuffle',100); % # of shuffle test in continuous process case
addParameter(p,'Fs',1000); % Only effective in time series data (Hz)
addParameter(p,'Lossfun','L2'); % Definition of loss function
parse(p,varargin{:});

if rank([ones(size(P,1),1), P]) <= size(P,2)
    error('The spatial information contains less dimensions than it shows.')
end

if iscell(Input)
    %% Point process case (multiunit, negative peak, maximal descent)
    n_channel = numel(Input);
    P = mat2cell(P,ones(n_channel,1),size(P,2));
    for i = 1:n_channel
        Input{i} = Input{i}(:); % make sure the alignment is correct
        P{i} = repmat(P{i}(:)',[numel(Input{i}),1]);
    end
    P = cell2mat(P);
    P1 = P(:,1);
    P2 = P(:,2);
    X = P;
    Y = cell2mat(Input(:));      
    % Check sample size
    [Xm,Xn] = size(X);
    if Xm<=Xn
        warning('Not enough event sample (under-determined system), not enough data to perform estimation.')
        V = nan(1,Xn);
        pValue = nan;
        return;
    end
    
    % Linear regression - by loss function
    switch p.Results.Lossfun
        case 'L2'
            XC = [X,ones(Xm,1)];            
            Beta = lscov(XC,Y); 
            Fstat = ((var(Y,0)-var(Y-XC*Beta,0))/(Xn)) ...
                   /(var(Y-XC*Beta,0)/(Xm-Xn-1));
            F_CDF = @(f,df1,df2) betainc(f/(f+df2/df1),df1/2,df2/2);
            try
                pValue = 1- F_CDF(Fstat,Xn,Xm-Xn-1); % F-test
            catch
                pValue = 1;
            end
        case 'L1'
            [Beta, E, E_s, pValue] = L1LinearRegression(X,Y); % Nested function
    end
elseif ismatrix(Input) && isnumeric(Input)
    %% Time series case
    W = Input;
    [n_channel, n_t] = size(W);
    T_vector = (1:n_t)/p.Results.Fs;
    P1 = P(:,1);
    P2 = P(:,2);
    P1 = repmat(P1,[1,n_t]);
    P2 = repmat(P2,[1,n_t]);
    P1 = P1(:);
    P2 = P2(:);
    % Prepare design matrix    
    X = [P1, P2]; 
    % Prepare response matrix % weight
    Y = repmat(T_vector,[n_channel,1]);
    Y = Y(:);
    W = W(:);
    % Check sample size
    [Xm,Xn] = size(X);
    if Xm<Xn
        warning('Not enough event sample (under-determined system), not enough data to perform estimation.')
        V = nan(1,Xn);
        return;
    end    
    switch p.Results.Lossfun
        case 'L2'
            [Beta, E, E_s, pValue] = L2LinearRegression(X, Y, W);            
        case 'L1'
            [Beta, E, E_s, pValue] = L1LinearRegression(X, Y, W);
    end
end
V = pinv(Beta(1:2));

%% Plotting the results
if p.Results.switch_plot
    % Plot the regression plane
    figure;
	subplot(2,2,1)
    if iscell(Input)
        scatter3(P1,P2,Y,'filled');hold on;        
    elseif ismatrix(Input) && isnumeric(Input)
        scatter3(P1,P2,Y,100*W/max(W(:)),'filled');hold on;
    end
    [P1U,P2U] = meshgrid(sort(unique(P1)),sort(unique(P2)));
    f = scatteredInterpolant(P1,P2,X*Beta(1:2) + Beta(end));
    Z = f(P1U,P2U);
    mesh(P1U,P2U,Z) %interpolated

%    legend('Data','Regression fit');
    xlabel('cm');ylabel('cm');zlabel('Second');

    % Plot the projection along the velocity axis
    subplot(2,2,3)
    P_v_axis = X*Beta(1:2)/norm(Beta(1:2));
    if iscell(Input)
        plot(P_v_axis,Y,'.');
    elseif ismatrix(Input) && isnumeric(Input)
        f = scatteredInterpolant(P_v_axis,Y,W);    
        P_v_axis_lin = linspace(min(P_v_axis),max(P_v_axis),n_channel);
        [P1U,P2U] = meshgrid(P_v_axis_lin,Y);    
        Z = f(P1U,P2U);
        imagesc(P_v_axis_lin,Y,Z);
    end
    hold on;
    plot(P_v_axis,[P1,P2]*Beta(1:2)+Beta(end));
    title('Projection along the velocity vector');
    xlabel('cm');
    ylabel('Second');colormap('hot')

    % Plot shuffled residual error distribution if L1-regression is used
    if strcmpi(p.Results.Lossfun,'L1')
        subplot(2,2,4)
        hist(E_s);hold on;
        plot([E E],get(gca,'YLim'),'r');hold off
        title(['Mean residual error, p = ' num2str(pValue)]);
    end    
end 


%% Nested function - L2 regression

    function [B, E, E_s, L2_p_shuffle] = L2LinearRegression(X,Y,Weight) 
        if nargin < 3
            Weight = 1;
        end
        % Determine size of predictor data 
        [n,~] = size(X); 
        % Add constant term
        X = [X,ones(n,1)];
        [B,~,E] = lscov(X,Y,Weight);
        % pValue
        if p.Results.n_shuffle
            % Shuffle the channel's spatial information
            E_s = nan(p.Results.n_shuffle,1);
                for i_shuffle = 1:p.Results.n_shuffle                 
                    [~,~,E_s(i_shuffle)] = lscov(Shuffle_Position(X),Y,Weight);
                end
                L2_p_shuffle = sum(E_s<=E)/p.Results.n_shuffle;
        else
            E_s = [];
            L2_p_shuffle = nan;          
        end  
    end

%% Nested function - L1 regression

    function [B, E, E_s, L1_p_shuffle] = L1LinearRegression(X,Y,Weight) 
        if nargin < 3
            Weight = 1;
        end
        % Determine size of predictor data 
        [n,m] = size(X); 
        % Add constant term
        X = [X,ones(n,1)];
        [B,E] = L1Optimize(X,Y,Weight);
        
        if p.Results.n_shuffle
            % Shuffle the channel's spatial information
            E_s = nan(p.Results.n_shuffle,1);
                for i_shuffle = 1:p.Results.n_shuffle                 
                    [~,E_s(i_shuffle)] = L1Optimize(Shuffle_Position(X),Y,Weight);
                end
                L1_p_shuffle = sum(E_s<=E)/p.Results.n_shuffle;
        else
            E_s = [];
            L1_p_shuffle = nan;
        end        
        function [Bout,Eout] = L1Optimize(Xin,Yin,Win)
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
                Bout = (repmat(w,[1 m+1]) .* Xin) \ (w .* Yin); 
                iter = iter + 1;
                if iter > 30;break;end
            end 
            % Report mean absolute deviation
            Eout = sum(w.* abs(Xin * BOld - Yin));
        end
    end

    function Xout = Shuffle_Position(Xin)
        [X_u,~,IB_X] = unique(Xin,'rows');
        Per = randperm(size(X_u,1));
        X_u_per = X_u(Per,:);
        Xout = X_u_per(IB_X,:);                          
    end

end
