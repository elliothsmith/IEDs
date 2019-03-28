% [Beta, SIGMA, pValue] = SpatialLinearRegression(response,P,varargin)
%
% 2-Dimensional linear regression to find whether spatial information is
% correlated with the observed response.  For detailed explanation, please
% refer to the 'Spatial linear regression methodology.docx'
%
% Input: 1) response: 
%           ** Important **
%           Input type can be cell or matrix, depends on whether they are
%           point process or time series data.  
%           [Data type: cell] 
%           For point process data.  One channel can have multiple events.
%           Unit: second
%           Commonly used data are:
%           - invasion time (each channel's seizure/burst invasion time)
%           - burst information, such as f_peak, peak_time ... etc.
%           - timestamps of unit activity
%           [Data type: matrix]
%           For time series data.  It should be a n_channel by n_t matrix.
%           This will launch weighted regression and shuffle test
%           validation.
%        2) P:
%           The position information. Format: an n_channel by n_dimension
%           matrix.  Currently I only allow n_dimension <= 2.
%        3) [Optional] parameter_name - parameter
%           alpha: significance level, default 0.05
%           switch_plot: (1/0) - plot the result or not
%           n_shuffle: number of shuffle test for time series data
%           Fs: sampling frequency for time series data - default:1000 Hz
%           Lossfun: (L1/L2), define the loss function in regression.
%                    Default: L2
%
% Output: V: Estimated velocity
%         H: Significance (1/0) given alpha
%            For point process data: F-test
%            For time series data: Shuffle test for residual error, notice
%            this is actually equivalent to F-test in shuffled data
%         pValue: p Value of F-test or shuffle test
%
% Detailed algorithm: Please see the word file.
%
% Called by SpatialAnalysis.m 
%           SpatialAnalysisLFP.m 
%           SpatialAnalysisInterictal.m 
%
% Author: Jyun-you Liou
% Final update: 2016/09/14 (This doesn't require statsitics toolbox)
% Now also report the estimator's covariance matrix 
% 
% [20181019EHS] TODO:: output df and test statistics
function [Beta, SIGMA, pValue] = SpatialLinearRegression(response,P,varargin)
p = inputParser;
addParameter(p,'alpha',0.05); % Significance level
addParameter(p,'switch_plot',0); 
   % In point process, it will plot data average at its physical position.
   % In continuous process, it will return shuffle test results.
addParameter(p,'n_shuffle',100); % # of shuffle test in continuous process case
addParameter(p,'Fs',1000); % Only effective in time series data (Hz)
addParameter(p,'Lossfun','L2'); % Definition of loss function
parse(p,varargin{:});

if iscell(response)
    %% Point process case
    n_channel = numel(response);
    P = mat2cell(P,ones(n_channel,1),size(P,2));
    for i = 1:n_channel
        response{i} = response{i}(:); % make sure the alignment is correct
        P{i} = repmat(P{i}(:)',[numel(response{i}),1]);
    end
    P = cell2mat(P);
    X = P;
    P1 = P(:,1);
    P2 = P(:,2);
    Y = cell2mat(response(:));  
    % PCA to check the dimensionality of P1, P2
    [C,S,L] = pca(X);
    if any((L/sum(L))<0.001)
        warning('Available channels are less than 2-dimensions.')
        Dim_singular = all(S < 10^-10); % Find the singular dimension
        X = S(:,~Dim_singular);
    end
    % Check sample size
    [Xm,Xn] = size(X);
    if Xm<=Xn
        warning('Under-determined system, not enough data to estimate.')
        V = nan(1,Xn);
        H = false;
        pValue = nan;
        return;
    end
    % Linear regression - by loss function
    switch p.Results.Lossfun
        case 'L2'
            XC = [ones(Xm,1),X];            
            [Beta,~,~,SIGMA] = lscov(XC,Y); 
            Fstat = ((var(Y,0)-var(Y-XC*Beta,0))/(Xn)) ...
                   /(var(Y-XC*Beta,0)/(Xm-Xn-1));
            F_CDF = @(f,df1,df2) betainc(f/(f+df2/df1),df1/2,df2/2);
            try
                pValue = 1- F_CDF(Fstat,Xn,Xm-Xn-1); % F-test
            catch
                pValue = 1;
            end
        case 'L1'
            [Beta, E, E_s, pValue, SIGMA] = L1LinearRegression(X,Y); % Nested function
    end
    H = pValue < p.Results.alpha;
    % Report beta
    if any((L/sum(L))<0.001) % For situation if there is spatial degeneracy
        B = Beta(2:end);
        B_rotated = zeros(numel(Dim_singular),1);
        B_rotated(~Dim_singular) = B;
        B = C*B_rotated;
    else
        B = Beta(2:end);
    end
elseif ismatrix(response) && isnumeric(response)
    %% Time series case
    W = response;
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
    % PCA to check the dimensionality of P1, P2
    [C,S,L] = pca(X);
    if any((L/sum(L))<0.001)
        warning('Available channels are less than 2-dimensions.')
        Dim_singular = all(S < 10^-10); % Find the singular dimension
        X = S(:,~Dim_singular);
    end    
    % Prepare response matrix % weight
    Y = repmat(T_vector,[n_channel,1]);
    Y = Y(:);
    W = W(:);
    % Check sample size
    [Xm,Xn] = size(X);
    if Xm<Xn
        warning('Under-determined system, not enough data to estimate.')
        V = nan(1,Xn);
        H = false;
        return;
    end    
    switch p.Results.Lossfun
        case 'L2'
            [Beta, E, E_s, pValue, SIGMA] = L2LinearRegression(X, Y, W);            
        case 'L1'
            [Beta, E, E_s, pValue, SIGMA] = L1LinearRegression(X, Y, W);
    end
    H = pValue <= p.Results.alpha;
    % Report beta
    if any((L/sum(L))<0.001) % For situation if there is spatial degeneracy
        B = Beta(2:end);
        B_rotated = zeros(numel(Dim_singular),1);
        B_rotated(~Dim_singular) = B;
        B = C*B_rotated;
    else
        B = Beta(2:end);
    end  
end
V = pinv(B);

%% Plotting the results
if p.Results.switch_plot
    % Plot the regression plane
    figure(gcf);
    subplot(2,2,[1 2]);
    if iscell(response)
        scatter3(P1,P2,Y,'filled');hold on;        
    elseif ismatrix(response) && isnumeric(response)
        scatter3(P1,P2,Y,100*W/max(W(:)),'filled');hold on;
    end
    if any((L/sum(L))<0.001) % 1D situation
        plot3(P1,P2,X*Beta(2:end) + Beta(1));
    else
        [P1U,P2U] = meshgrid(sort(unique(P1)),sort(unique(P2)));
        f = scatteredInterpolant(P1,P2,X*Beta(2:end) + Beta(1));
        Z = f(P1U,P2U);
        mesh(P1U,P2U,Z) %interpolated
    end
    legend('Data','Regression fit');
    xlabel('cm');ylabel('cm');zlabel('Second');
    
    % Plot the projection along the velocity axis
    subplot(2,2,3)
    P_v_axis = X*Beta(2:end)/norm(Beta(2:end));
    if iscell(response)
        plot(P_v_axis,Y,'.');
    elseif ismatrix(response) && isnumeric(response)
        f = scatteredInterpolant(P_v_axis,Y,W);    
        P_v_axis_lin = linspace(min(P_v_axis),max(P_v_axis),n_channel);
        [P1U,P2U] = meshgrid(P_v_axis_lin,Y);    
        Z = f(P1U,P2U);
        imagesc(P_v_axis_lin,Y,Z);
    end
    hold on;
    plot(P_v_axis,[P1,P2]*Beta(2:end)+Beta(1));
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
    %waitforbuttonpress;
end 


%% Nested function - L2 regression

    function [B, E, E_s, L2_p_shuffle, SIGMA] = L2LinearRegression(X,Y,Weight) 
        if nargin < 3
            Weight = 1;
        end
        % Determine size of predictor data 
        [n,~] = size(X); 
        % Add constant term
        X = [ones(n,1),X];
        [B,~,E,SIGMA] = lscov(X,Y,Weight);
        if p.Results.n_shuffle
            % Shuffle the channel's spatial information
            E_s = nan(p.Results.n_shuffle,1);
                for i_shuffle = 1:p.Results.n_shuffle                 
                    [~,~,E_s(i_shuffle)] = lscov(Shuffle_Position(X),Y,Weight);
                end
                L2_p_shuffle = sum(E_s<=E)/p.Results.n_shuffle;
        else
            L2_p_shuffle = nan;          
        end     
    end

%% Nested function - L1 regression

    function [B, E, E_s, L1_p_shuffle, SIGMA] = L1LinearRegression(X,Y,Weight) 
        if nargin < 3
            Weight = 1;
        end
        % Determine size of predictor data 
        [n,m] = size(X); 
        % Add constant term
        X = [ones(n,1),X];
        [B,E,SIGMA] = L1Optimize(X,Y,Weight);
        
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
        function [Bout,Eout,SIGMA] = L1Optimize(Xin,Yin,Win)
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
            % Calculate SIGMA
            [~,~,~,SIGMA] = lscov(Xin,Yin,w);
        end
    end

    function Xout = Shuffle_Position(Xin)
        [X_u,~,IB_X] = unique(Xin,'rows');
        Per = randperm(size(X_u,1));
        X_u_per = X_u(Per,:);
        Xout = X_u_per(IB_X,:);                          
    end

end