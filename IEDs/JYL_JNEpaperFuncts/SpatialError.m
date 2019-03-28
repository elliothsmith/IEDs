% This help you to calculate the error of predictors made by either gamma
% or LFP method
%
% Input: 1. time & estimator from spikes (data1) - as golden standard
%        2. time & estimator from another method (e.g. gamma/LFP) (data2)
%           allowed for multiple comparison in a structure
%
% Output: Error in direction and speed
%         varargout: McNemar test for individual 2 by 2 table
%
% This function calls circular_stat toolbox
%
% Author: Jyun-you Liou
% Final update: 2015/11/25

function [ S_err, D_err, H_err, pMTest, Sen, Spe ] = SpatialError( data1, data2, varargin )

p = inputParser;
addParameter(p,'dT',0.05); % The error of timing of a single burst detected by the two different methods
parse(p,varargin{:});
dT = p.Results.dT; 
n_method = numel(data2);
% From data1 (spike data)
T1 = data1.T(:); % timing of the burst
H1 = data1.Significance(:); % Significance
S1 = sqrt(sum((data1.V).^2,2)); % Speed, unit: cm/second
D1 = data1.direction(:); % Direction in radian
C1 = false(numel(T1),numel(data2)); % whether you find corresponding episodes or not
TC = nan(numel(T1),numel(data2)); % The corresponding episode timing
HC = false(numel(T1),numel(data2)); % The corresponding episode significance
SC = nan(numel(T1),numel(data2)); % The corresponding episode speed
DC = nan(numel(T1),numel(data2)); % The corresponding episode direction

for i= 1:numel(T1);
    for m = 1:numel(data2)
        % Find the most nearby burst
        [T_error,LOC] = min(abs(data2(m).T-T1(i)));
        if T_error < dT
            C1(i,m) = true;
            TC(i,m) = data2(m).T(LOC);
            HC(i,m) = data2(m).Significance(LOC);          
            SC(i,m) = norm(data2(m).V(LOC,:));
            DC(i,m) = data2(m).direction(LOC);
        end
    end
end

% Evaluate performance
H_success = logical(bsxfun(@times,HC,H1));
H_sen = logical(bsxfun(@times,HC,H1));
Sen = sum(H_sen) ./ sum(H1);
H_spe = logical(bsxfun(@times,~HC,~H1));
Spe = sum(H_spe) ./ sum(~H1);
H_err = sum(~bsxfun(@xor,H1,HC)) ./ size(H1,1);

S_err = bsxfun(@minus,SC,S1);
D_err(:,:,1) = bsxfun(@minus,DC,D1);
D_err(:,:,2) = bsxfun(@minus,DC,D1) + 2*pi;
D_err(:,:,3) = bsxfun(@minus,DC,D1) - 2*pi;
Selector = D_err<pi & D_err>-pi;
D_err = D_err .* Selector;
D_err = squeeze(sum(D_err,3));



figure(1)
subplot(2,1,1);
plot(T1,S_err,'.');hold on;
T_plot = repmat(T1,[1 n_method]);
plot(T_plot(H_success),S_err(H_success),'rx');
title('Error in speed');
% for m = 1:n_method
%     text{m} = ['Method ' data2(m).method ': ' num2str(round(1000*H_err(m))/10) ' % agreement'];
% end
% text{n_method +1}='Detected discharges.';
% legend(text,'location','best');
% subplot(2,1,2);
% plot(T1,D_err,'.');hold on;
% plot(T_plot(H_success),D_err(H_success),'rx');
% title('Error in direction')
% legend(text,'location','best');
% clear text;

S_err_hist = S_err;
D_err_hist = D_err;
S_err_hist(~H_success) = nan;
D_err_hist(~H_success) = nan;
% We don't use nanstd nanmean here so the user doeasn't need the financial toolbox
S_err_mean = sum (S_err .* H_success) ./ sum(H_success); 

for m = 1:n_method
    Temp = S_err_hist(:,m);
    Temp(isnan(Temp)) = [];
    S_err_std(m) = std(Temp);
    S_err_median(m) = median(S_err(H_success(:,m),m));
    S_err_quantile(m) = quantile(S_err(H_success(:,m),m),0.75)-quantile(S_err(H_success(:,m),m),0.25);    
end

figure(2)
subplot(2,1,1)
hist(S_err_hist,-50:3:100);
ylabel('# of discharges');
xlabel('cm/sec');
title('Error in speed');
% for m = 1:n_method
%     text{m} = ['Method ' data2(m).method ': ' num2str(round(1000*H_err(m))/10) ' % agreement' ...
%                sprintf('\nError median = ') num2str(S_err_median(m)) ...
%                sprintf('\nError quantile = ') num2str(S_err_quantile(m))];
% end
% legend(text,'location','best');
% clear text;

subplot(2,1,2)
hist(D_err_hist,40);
title('Error in direction');
ylabel('# of discharges')
xlabel('radian');
xlim([-pi pi]);
for m = 1:n_method
    Temp = D_err_hist(:,m);
    Temp(isnan(Temp)) = [];
    D_err_mean(m) = circ_mean(Temp);
    D_err_std(m) = circ_std(Temp);
end

% for m = 1:numel(data2)
%     text{m} = ['Method ' data2(m).method ': ' num2str(round(1000*H_err(m))/10) ' % agreement' ...
%                sprintf('\nError mean = ') num2str(D_err_mean(m)) ...
%                sprintf('\nError std = ') num2str(D_err_std(m))];
% end
% legend(text,'location','best');

if nargout > 3
    pMTest = zeros(numel(data2));
    for i = 1:numel(data2)
        for j = 1:i-1
            [~,pMTest(i,j)] = testcholdout(HC(:,i),HC(:,j),H1);
        end
    end
    pMTest = pMTest + pMTest';
end

end