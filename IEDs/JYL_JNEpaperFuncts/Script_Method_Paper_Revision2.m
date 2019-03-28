%% For revision - regenerate all data
Fs = 1000; % Hz
[b1,a1] = fir1(90,50/(Fs/2));
[b2,a2] = fir1(90,[80 150]/(Fs/2));
n_ch = [(3:8).^2,inf];
load('RevisionData.mat', 'Record_combined')
Trec = Record_combined(2,1).T;

%% Methods - neg-peak, max-descent, gamma, xcorr
Beta = nan(3,numel(Trec),4,2,numel(n_ch)); % Beta = 3 by n_burst by n_method by L2/L1 by 7 by condition
Sigma = nan(3,3,numel(Trec),4,2,numel(n_ch)); % Sigma = 3 by n_burst by n_method by L2/L1 by 7 by condition
pValue = nan(numel(Trec),4,2,numel(n_ch)); % pValue = n_burst by n_method by L2/L1 by condition
for iter = 7
    n_channel = n_ch(iter);
    load('C5.mat'); % The T_burst is correct, as C5_new
    clear P;
    
    % map
    ch = sort(map(map>0));
    for i = numel(ch):-1:1
        [P(i,1),P(i,2)] = find(map == ch(i));
    end
    P = P*0.04; % cm

    % LFP
    LFP = LFP(ch,:);
    LFP_low = filtfilt(b1,a1,LFP')';
    LFP_gamma = filtfilt(b2,a2,LFP')';
    LFP_power = (abs(hilbert(LFP_gamma')')).^2;
    
    % Run function
    for i = 1:numel(T_burst)
        tic;
        p = round(T_burst(i) * 1000);
        LFP_local_low = LFP_low(:,p + (-100:100));
        LFP_local_power = LFP_power(:,p + (-30:30));
        P_local = P;
        if n_channel < inf
            ch_selection = randperm(size(LFP_local_low,1));
            ch_selection = ch_selection(1:n_channel);
            LFP_local_low = LFP_local_low(ch_selection,:);
            LFP_local_power = LFP_local_power(ch_selection,:);            
            P_local = P_local(ch_selection,:);
        end
        
        % Neg-peak n_burst by n_method by L2/L1 by 7 by condition
        [~,p_neg] = min(LFP_local_low(:,51:151),[],2);
        t_neg = p_neg / Fs;t_neg = mat2cell(t_neg,ones(numel(t_neg),1));       
        [Beta(:,i,1,1,iter), Sigma(:,:,i,1,1,iter), pValue(i,1,1,iter)] = SpatialLinearRegression( t_neg, P_local );
        [Beta(:,i,1,2,iter), Sigma(:,:,i,1,2,iter), pValue(i,1,2,iter)] = SpatialLinearRegression( t_neg, P_local, 'Lossfun', 'L1' );

        % Max-descent
        [~,p_max_descent] = min(diff(LFP_local_low(:,51:151),1,2),[],2);
        t_max_descent = p_max_descent / Fs;t_max_descent = mat2cell(t_max_descent,ones(numel(t_max_descent),1));
        [Beta(:,i,2,1,iter), Sigma(:,:,i,2,1,iter), pValue(i,2,1,iter)] = SpatialLinearRegression( t_max_descent, P_local );
        [Beta(:,i,2,2,iter), Sigma(:,:,i,2,2,iter), pValue(i,2,2,iter)] = SpatialLinearRegression( t_max_descent, P_local, 'Lossfun', 'L1' );

        % Gamma
        [Beta(:,i,3,1,iter), Sigma(:,:,i,3,1,iter), pValue(i,3,1,iter)] = SpatialLinearRegression( LFP_local_power, P_local );
        [Beta(:,i,3,2,iter), Sigma(:,:,i,3,2,iter), pValue(i,3,2,iter)] = SpatialLinearRegression( LFP_local_power, P_local, 'Lossfun', 'L1' );
        
        % Xcorr n_burst by n_method by L2/L1 by 7 by condition
        [Beta(2:3,i,4,1,iter), Sigma(2:3,2:3,i,4,1,iter), pValue(i,4,1,iter)] = SpatialLinearRegressionXcorr( LFP_local_low, P_local );
        [Beta(2:3,i,4,2,iter), Sigma(2:3,2:3,i,4,2,iter), pValue(i,4,2,iter)] = SpatialLinearRegressionXcorr( LFP_local_low, P_local, 'Lossfun', 'L1' );

        toc;
    end
    %% Nu
    load('Nu.mat');
    clear P;

    % map
    ch = sort(map(map>0));
    for i = numel(ch):-1:1
        [P(i,1),P(i,2)] = find(map == ch(i));
    end
    P = P*0.04;

    % LFP
    LFP = LFP(ch,:);
    LFP_low = filtfilt(b1,a1,LFP')';
    LFP_gamma = filtfilt(b2,a2,LFP')';
    LFP_power = (abs(hilbert(LFP_gamma')')).^2;
    
    % Run function
    for i = 1:numel(T_burst)
        tic;
        p = round(T_burst(i) * 1000);
        LFP_local_low = LFP_low(:,p + (-100:100));
        LFP_local_power = LFP_power(:,p + (-30:30));        
        P_local = P;
        if n_channel < inf
            ch_selection = randperm(size(LFP_local_low,1));
            ch_selection = ch_selection(1:n_channel);
            LFP_local_low = LFP_local_low(ch_selection,:);
            LFP_local_power = LFP_local_power(ch_selection,:);            
            P_local = P_local(ch_selection,:);
        end
        
        % Neg-peak n_burst by n_method by L2/L1 by 7 by condition
        [~,p_neg] = min(LFP_local_low(:,51:151),[],2);
        t_neg = p_neg / Fs;t_neg = mat2cell(t_neg,ones(numel(t_neg),1));       
        [Beta(:,i+1029,1,1,iter), Sigma(:,:,i+1029,1,1,iter), pValue(i+1029,1,1,iter)] = SpatialLinearRegression( t_neg, P_local );
        [Beta(:,i+1029,1,2,iter), Sigma(:,:,i+1029,1,2,iter), pValue(i+1029,1,2,iter)] = SpatialLinearRegression( t_neg, P_local, 'Lossfun', 'L1' );

        % Max-descent
        [~,p_max_descent] = min(diff(LFP_local_low(:,51:151),1,2),[],2);
        t_max_descent = p_max_descent / Fs;t_max_descent = mat2cell(t_max_descent,ones(numel(t_max_descent),1));
        [Beta(:,i+1029,2,1,iter), Sigma(:,:,i+1029,2,1,iter), pValue(i+1029,2,1,iter)] = SpatialLinearRegression( t_max_descent, P_local );
        [Beta(:,i+1029,2,2,iter), Sigma(:,:,i+1029,2,2,iter), pValue(i+1029,2,2,iter)] = SpatialLinearRegression( t_max_descent, P_local, 'Lossfun', 'L1' );

        % Gamma
        [Beta(:,i+1029,3,1,iter), Sigma(:,:,i+1029,3,1,iter), pValue(i+1029,3,1,iter)] = SpatialLinearRegression( LFP_local_power, P_local );
        [Beta(:,i+1029,3,2,iter), Sigma(:,:,i+1029,3,2,iter), pValue(i+1029,3,2,iter)] = SpatialLinearRegression( LFP_local_power, P_local, 'Lossfun', 'L1' );
        
        % Xcorr n_burst by n_method by L2/L1 by 7 by condition
        [Beta(2:3,i+1029,4,1,iter), Sigma(2:3,2:3,i+1029,4,1,iter), pValue(i+1029,4,1,iter)] = SpatialLinearRegressionXcorr( LFP_local_low, P_local );
        [Beta(2:3,i+1029,4,2,iter), Sigma(2:3,2:3,i+1029,4,2,iter), pValue(i+1029,4,2,iter)] = SpatialLinearRegressionXcorr( LFP_local_low, P_local, 'Lossfun', 'L1' );
        toc
    end

    %% n_condition by n_method*L2/L1
    str = {'Neg-L2','MaxDesc-L2','Gamma-L2','Xcorr-L2','Neg-L1','MaxDesc-L1','Gamma-L1','Xcorr-L1'};
    for k = 1:8
        Rec(iter,k).T = Trec;
        Rec(iter,k).B   = Beta (2:3,  :,   mod(k-1,4)+1, ceil(k/4), iter);        
        Rec(iter,k).Cov = Sigma(2:3,2:3,:, mod(k-1,4)+1, ceil(k/4), iter);
        Rec(iter,k).p =  pValue(        :, mod(k-1,4)+1, ceil(k/4), iter);
        for l = numel(Rec(iter,k).p):-1:1
            V(l,:) = pinv(Rec(iter,k).B(:,l));
        end
        Rec(iter,k).V = V;
        Rec(iter,k).direction = angle(V(:,1) + sqrt(-1)*V(:,2));        
        Rec(iter,k).Significance = Rec(iter,k).p<0.05;
        Rec(iter,k).method = str{k};        
    end

end
save('NewResults.mat','Rec','-v7.3')

%% Produce comparison 
% load('D:\MATLAB\Eilepsy Codes\Data\method paper data\RevisionData.mat', 'Record_combined')
%%
% load('D:\MATLAB\Eilepsy Codes\Data\method paper data\RevisionData.mat', 'Record_LE_combined_all')
% [Record_LE_combined_all.p] = Record_LE_combined_all.pValue;
% Record_LE_combined_all = rmfield(Record_LE_combined_all,'pValue');
% RecAll = [Record_LE_combined_all,RecXcorr,RecAvg];
% 
%%
% load('RevisionAllResults.mat');
% str = {'L2Neg', 'L2MaxD', 'L2Envelope', 'L2SimpleAvg', 'L1Neg', 'L1MaxD', ...
%        'L1Envelope', 'L1SimpleAvg', 'Xcorr', 'XcorrWeighted', 'L2WeightedAvg', ...
%        'L2OffsetAvg', 'L1WeightedAvg', 'L1OffsetAvg'};
% 
% [ DD, DS, DH, pD, pS, pH, dD, dS, DR ] = CompareMethods( RecAll, Record_MUA_combined(1,1), [0.25 0.5 0.75]  );
% %%
% f = figure;A = axes;
% % L = plot(dS(:,:,2));
% L = plot(DR);
% 
% for i = 1:numel(L)
%     L(i).Color = rand(3,1);
% end
% legend(str)
