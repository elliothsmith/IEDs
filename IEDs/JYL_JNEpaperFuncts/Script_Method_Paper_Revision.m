%% For revision - regenerate all data
Fs = 1000; % Hz
[b1,a1] = fir1(90,50/(Fs/2));
[b2,a2] = fir1(90,[80 150]/(Fs/2));
n_ch = [(3:8).^2,inf];
load('D:\MATLAB\Eilepsy Codes\Data\method paper data\RevisionData.mat', 'Record_combined')
Trec = Record_combined(2,1).T;

%% Methods - neg-peak, max-descent, gamma, xcorr
Beta = nan(2,numel(Trec),4,2,numel(n_ch)); % Beta = 2 by n_burst by n_method by L2/L1 by 7 by condition
Sigma = nan(2,2,numel(Trec),4,2,numel(n_ch)); % Sigma = 2 by n_burst by n_method by L2/L1 by 7 by condition
pValue = nan(numel(Trec),4,2,numel(n_ch)); % pValue = n_burst by n_method by L2/L1 by condition
for iter = 1:7
    n_channel = n_ch(iter);
    load('D:\MATLAB\Eilepsy Codes\Data\method paper data\C5.mat'); % The T_burst is correct, as C5_new

    % map
    ch = sort(map(map>0));
    for i = numel(ch):-1:1
        [P(i,1),P(i,2)] = find(map == ch(i));
    end
    P = P*0.04; % cm

    % LFP
    LFP = LFP(ch,:);
    LFP_low = filtfilt(b1,a1,LFP')';
    LFP_gamma = filtfilt(b1,a1,LFP')';
    LFP_power = (hilbert(LFP_gamma')').^2;
    
    % Run function
    for i = 1:numel(T_burst)
        tic;
        p = round(T_burst(i) * 1000);
        LFP_local_low = LFP_low(:,p + (-100:100));
        LFP_local_gamma = LFP_gamma(:,p + (-50:50));        
        P_local = P;
        if n_channel < inf
            ch_selection = randperm(size(LFP_local_low,1));
            ch_selection = ch_selection(1:n_channel);
            LFP_local_low = LFP_local_low(ch_selection,:);
            LFP_power = LFP_power(ch_selection,:);            
            P_local = P_local(ch_selection,:);
        end
        
        % Neg-peak n_burst by n_method by L2/L1 by 7 by condition
        [Beta(:,i,1,1,iter), Sigma(:,:,i,1,1,iter), pValue(i,1,1,iter)] = SpatialLinearRegression( LFP_local_low, P_local );
        [Beta(:,i,1,2,iter), Sigma(:,:,i,1,2,iter), pValue(i,1,2,iter)] = SpatialLinearRegression( LFP_local_low, P_local, 'Lossfun', 'L1' );

        % Max-descent
        [Beta(:,i,2,1,iter), Sigma(:,:,i,2,1,iter), pValue(i,2,1,iter)] = SpatialLinearRegression( LFP_local_low, P_local );
        [Beta(:,i,2,2,iter), Sigma(:,:,i,2,2,iter), pValue(i,2,2,iter)] = SpatialLinearRegression( LFP_local_low, P_local, 'Lossfun', 'L1' );

        

        if strcmpi(Method,'Xcorr')
            [ Beta(:,i), pValue(i) ] = SpatialLinearRegressionXcorr( LFP_local_low, P_local );
            [ Beta1(:,i), pValue1(i) ] = SpatialLinearRegressionXcorr( LFP_local_low, P_local, 'Lossfun', 'L1' );
        elseif strcmpi(Method,'averaged')
            LFP_local_low = LFP_local_low(:,51:150);            
            [~,p_neg] = min(LFP_local_low,[],2);
            [~,p_max_descent] = min(diff(LFP_local_low,1,2),[],2);
            t_neg = p_neg / Fs;
            t_max_descent = p_max_descent / Fs;
            % L1
            [V_n,~,pValue_n,Sigma_n] = SpatialLinearRegression(mat2cell(t_neg,ones(numel(t_neg),1)),P_local);
            [V_m,~,pValue_m,Sigma_m] = SpatialLinearRegression(mat2cell(t_max_descent,ones(numel(t_max_descent),1)),P_local);
            Beta_n = pinv(V_n);
            Beta_m = pinv(V_m);
            Beta(:,i) = (inv(Sigma_n(1:2,1:2)) + inv(Sigma_m(1:2,1:2))) \ (Sigma_n(1:2,1:2)\Beta_n + Sigma_m(1:2,1:2)\Beta_m);
            pValue(i) = min([pValue_n,pValue_m]);
            % L1
            [V_n,~,pValue_n,Sigma_n] = SpatialLinearRegression(mat2cell(t_neg,ones(numel(t_neg),1)),P_local, 'Lossfun', 'L1');
            [V_m,~,pValue_m,Sigma_m] = SpatialLinearRegression(mat2cell(t_max_descent,ones(numel(t_max_descent),1)),P_local, 'Lossfun', 'L1');
            Beta_n = pinv(V_n);
            Beta_m = pinv(V_m);
            Beta1(:,i) = (inv(Sigma_n(1:2,1:2)) + inv(Sigma_m(1:2,1:2))) \ (Sigma_n(1:2,1:2)\Beta_n + Sigma_m(1:2,1:2)\Beta_m);
            pValue1(i) = min([pValue_n,pValue_m]);
            
            % Another way to do it
            t_neg = t_neg - mean(t_neg);
            t_max_descent = t_max_descent - mean(t_max_descent);
            [V,~,pValue_a(i)] = SpatialLinearRegression(mat2cell([t_neg;t_max_descent],ones(numel(t_neg)+numel(t_max_descent),1)),[P_local;P_local]);
            Beta_a(:,+i) = pinv(V);
            [V,~,pValue_a1(i)] = SpatialLinearRegression(mat2cell([t_neg;t_max_descent],ones(numel(t_neg)+numel(t_max_descent),1)),[P_local;P_local], 'Lossfun', 'L1');
            Beta_a1(:,+i) = pinv(V);    
            
        end
        toc;
    end

    %% Nu
    load('D:\MATLAB\Eilepsy Codes\Data\method paper data\Nu.mat');

    % map
    ch = sort(map(map>0));
    for i = numel(ch):-1:1
        [P(i,1),P(i,2)] = find(map == ch(i));
    end
    P = P*0.04;

    % LFP
    LFP = LFP(ch,:);
    LFP = filtfilt(b,a,LFP')';

    % Run function
    for i = 1:numel(T_burst)
        tic;
        p = round(T_burst(i) * 1000);
        LFP_local_low = LFP(:,p + (-100:100));
        P_local = P;
        if n_channel < inf
            ch_selection = randperm(size(LFP_local_low,1));
            ch_selection = ch_selection(1:n_channel);
            LFP_local_low = LFP_local_low(ch_selection,:);
            P_local = P_local(ch_selection,:);           
        end
        if strcmpi(Method,'Xcorr')
            [ Beta(:,1029+i), pValue(1029+i) ] = SpatialLinearRegressionXcorr( LFP_local_low, P_local );
            [ Beta1(:,1029+i), pValue1(1029+i) ] = SpatialLinearRegressionXcorr( LFP_local_low, P_local, 'Lossfun', 'L1' );
        elseif strcmpi(Method,'averaged')
            LFP_local_low = LFP_local_low(:,51:150);
            [~,p_neg] = min(LFP_local_low,[],2);
            [~,p_max_descent] = min(diff(LFP_local_low,1,2),[],2);
            t_neg = p_neg / Fs;
            t_max_descent = p_max_descent / Fs;
            % L2
            [V_n,~,pValue_n,Sigma_n] = SpatialLinearRegression(mat2cell(t_neg,ones(numel(t_neg),1)),P_local);
            [V_m,~,pValue_m,Sigma_m] = SpatialLinearRegression(mat2cell(t_max_descent,ones(numel(t_max_descent),1)),P_local);
            Beta_n = pinv(V_n);
            Beta_m = pinv(V_m);
            Beta(:,1029+i) = (inv(Sigma_n(1:2,1:2)) + inv(Sigma_m(1:2,1:2))) \ (Sigma_n(1:2,1:2)\Beta_n + Sigma_m(1:2,1:2)\Beta_m);
            pValue(1029+i) = min([pValue_n,pValue_m]);
            % L1
            [V_n,~,pValue_n,Sigma_n] = SpatialLinearRegression(mat2cell(t_neg,ones(numel(t_neg),1)),P_local, 'Lossfun', 'L1');
            [V_m,~,pValue_m,Sigma_m] = SpatialLinearRegression(mat2cell(t_max_descent,ones(numel(t_max_descent),1)),P_local, 'Lossfun', 'L1');
            Beta_n = pinv(V_n);
            Beta_m = pinv(V_m);
            Beta1(:,1029+i) = (inv(Sigma_n(1:2,1:2)) + inv(Sigma_m(1:2,1:2))) \ (Sigma_n(1:2,1:2)\Beta_n + Sigma_m(1:2,1:2)\Beta_m);
            pValue1(1029+i) = min([pValue_n,pValue_m]);
            
            % Another way to do it
            t_neg = t_neg - mean(t_neg);
            t_max_descent = t_max_descent - mean(t_max_descent);
            [V,~,pValue_a(1029+i)] = SpatialLinearRegression(mat2cell([t_neg;t_max_descent],ones(numel(t_neg)+numel(t_max_descent),1)),[P_local;P_local]);
            Beta_a(:,1029+i) = pinv(V);
            [V,~,pValue_a1(1029+i)] = SpatialLinearRegression(mat2cell([t_neg;t_max_descent],ones(numel(t_neg)+numel(t_max_descent),1)),[P_local;P_local], 'Lossfun', 'L1');
            Beta_a1(:,1029+i) = pinv(V);            
        end
        toc;
    end


    %% 
    if strcmpi(Method,'xcorr') 
        Beta = -Beta;
        Beta1 = -Beta1;
    end
    
    for i = 1:numel(pValue)
        V(i,:) = pinv(Beta(:,i));
        V1(i,:) = pinv(Beta1(:,i));
        if strcmpi(Method,'averaged')
            V_a(i,:) = pinv(Beta_a(:,i));
            V_a1(i,:) = pinv(Beta_a1(:,i));            
        end
    end
    
    %%
    if strcmpi(Method,'xcorr')
        RecXcorr(iter,1).T = Trec;
        RecXcorr(iter,1).V = V;
        RecXcorr(iter,1).direction = angle(V(:,1) + sqrt(-1)*V(:,2));
        RecXcorr(iter,1).p = pValue';
        RecXcorr(iter,1).Significance = pValue'<0.05;
        RecXcorr(iter,1).B = Beta;
        RecXcorr(iter,1).method = 'xcorr';

        RecXcorr(iter,2).T = Trec;
        RecXcorr(iter,2).V = V1;
        RecXcorr(iter,2).direction = angle(V1(:,1) + sqrt(-1)*V1(:,2));
        RecXcorr(iter,2).p = pValue1';
        RecXcorr(iter,2).Significance = pValue1' < 0.05;
        RecXcorr(iter,2).B = Beta1;
        RecXcorr(iter,2).method = 'xcorr-L1';

    elseif strcmpi(Method,'averaged')
        RecAvg(iter,1).T = Trec;
        RecAvg(iter,1).V = V;
        RecAvg(iter,1).direction = angle(V(:,1) + sqrt(-1)*V(:,2));
        RecAvg(iter,1).p = pValue';
        RecAvg(iter,1).Significance = pValue'<0.025;
        RecAvg(iter,1).B = Beta;
        RecAvg(iter,1).method = 'avg';

        RecAvg(iter,2).T = Trec;
        RecAvg(iter,2).V = V1;
        RecAvg(iter,2).direction = angle(V1(:,1) + sqrt(-1)*V1(:,2));
        RecAvg(iter,2).p = pValue1';
        RecAvg(iter,2).Significance = pValue1' < 0.025;
        RecAvg(iter,2).B = Beta1;
        RecAvg(iter,2).method = 'avg-L1';

        RecAvg(iter,3).T = Trec;
        RecAvg(iter,3).V = V_a;
        RecAvg(iter,3).direction = angle(V_a(:,1) + sqrt(-1)*V_a(:,2));
        RecAvg(iter,3).p = pValue_a';
        RecAvg(iter,3).Significance = pValue_a'<0.025;
        RecAvg(iter,3).B = Beta_a;
        RecAvg(iter,3).method = 'avg-alt';

        RecAvg(iter,4).T = Trec;
        RecAvg(iter,4).V = V_a1;
        RecAvg(iter,4).direction = angle(V_a1(:,1) + sqrt(-1)*V_a1(:,2));
        RecAvg(iter,4).p = pValue_a1';
        RecAvg(iter,4).Significance = pValue_a1' < 0.025;
        RecAvg(iter,4).B = Beta_a1;
        RecAvg(iter,4).method = 'avg-alt-L1';        
        
    end

end

% save('XcorrResults.mat','RecXcorr','-v7.3')

%% Produce comparison 
load('D:\MATLAB\Eilepsy Codes\Data\method paper data\RevisionData.mat', 'Record_combined')
%%
% load('D:\MATLAB\Eilepsy Codes\Data\method paper data\RevisionData.mat', 'Record_LE_combined_all')
% [Record_LE_combined_all.p] = Record_LE_combined_all.pValue;
% Record_LE_combined_all = rmfield(Record_LE_combined_all,'pValue');
% RecAll = [Record_LE_combined_all,RecXcorr,RecAvg];
% 
%%
load('RevisionAllResults.mat');
str = {'L2Neg', 'L2MaxD', 'L2Envelope', 'L2SimpleAvg', 'L1Neg', 'L1MaxD', ...
       'L1Envelope', 'L1SimpleAvg', 'Xcorr', 'XcorrWeighted', 'L2WeightedAvg', ...
       'L2OffsetAvg', 'L1WeightedAvg', 'L1OffsetAvg'};

[ DD, DS, DH, pD, pS, pH, dD, dS, DR ] = CompareMethods( RecAll, Record_MUA_combined(1,1), [0.25 0.5 0.75]  );
%%
f = figure;A = axes;
% L = plot(dS(:,:,2));
L = plot(DR);

for i = 1:numel(L)
    L(i).Color = rand(3,1);
end
legend(str)
%% Try transport equation
% ut + cux = 0, ok give up

% Ut = diff(LFP_local,1,2); % uv/ms
% % For each channel, find its gradient
% for i = 1:size(P_local,1)
%     
% end
% Ux = bsxfun(@rdivide, diff(LFP_local,1), diff(P_local(:,1),1)); % uv/cm
% Uy = bsxfun(@rdivide, diff(LFP_local,1), diff(P_local(:,2),1));
% 
% ut = Ut(1:end-1,:);ut = ut(:);
% ux = Ux(:,1:end-1);ux = ux(:);
% uy = Uy(:,1:end-1);uy = uy(:);
% 
% ut(isinf(ut)) = nan;
% ux(isinf(ux)) = nan;
% uy(isinf(uy)) = nan;
% V = [ux, uy] \ ut;