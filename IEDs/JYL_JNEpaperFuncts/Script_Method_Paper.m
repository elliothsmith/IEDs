%% Fig 2
% load('/mnt/mfs/home/NIMASTER/jl4115/Epilepsy Codes/data/Human data/method paper data/C5.mat')
load('D:\MATLAB\Eilepsy Codes\Data\method paper data\C5.mat')
%% 2-1 & 4-1
c1 = [    0.4940    0.1840    0.5560];
c2 = [     0.4660    0.6740    0.1880];
c3 =  [        0    0.4470    0.7410];
c4 =  [    0.3010    0.7450    0.9330];
COLOR = [c1;c2;c3;c4];

Fs = 1000;
chan = 31;
ch = map(map>0);
%%
timestamp_to_raster(Ts(chan));
T_v = (1:size(LFP,2))/Fs;
[b,a]=butter(4,[70,200]/(Fs/2));
LFP_plot = filtfilt(b,a,LFP(chan,:))/100;
plot(T_v,LFP_plot);hold on;
% plot(T_v,abs(hilbert(LFP_plot)));hold on;
% legend({'MUA','High gamma (70~200 Hz) LFP','Instant amplitude of high gamma activity'},'location','best')
[~,neg_peak] = findpeaks(-LFP_plot,'minpeakheight',1);
p1 = plot(T_v(neg_peak),LFP_plot(neg_peak),'o');
p1.Color = c1;
[~,max_descent] = findpeaks(-diff(LFP_plot),'minpeakdistance',100);
p2 = plot(T_v(max_descent),LFP_plot(max_descent),'o');
p2.Color = c2;
legend({'MUA','LFP','Negative peak','Maximal descent'},'location','best');
ylabel('mV');xlabel('Second')

%% 2-2 4-2 Cross-correlation
% define all episodes
LFP_mean = mean(abs(hilbert(filtfilt(b,a,LFP'))),2);
% A_gamma = mean(abs(hilbert(filtfilt(b,a,LFP(ch,:)'))'));
% Ts_all = cell2mat(Ts(ch)');
% [~,~,CI] = CrossCorrelation({Ts_all},A_gamma,'nboot',200,'Fs',1000,'dT',0.05,'Fs_ts',2000);
[~,LOCS] = findpeaks(LFP_mean,'minpeakheight',LFP_threshold,'minpeakdistance',100);
Ts = cellfun(@(x) [x(x>15 & x<85);x(x>165 & x<250); x(x>315 & x<425)],Ts,'Uniformoutput',0);
chan = 31;
mask = 0*LFP(chan,:);mask(LOCS)=1;mask = conv(mask,ones(100,1),'same');
LFP_neg_peak = -LFP(chan,:) .* mask;
[~,T_neg]=findpeaks(LFP_neg_peak,'minpeakdistance',100);T_neg = T_neg/Fs;
dLFP_dt_minus = diff(LFP(chan,:));
dLFP_dt_minus = dLFP_dt_minus .* mask(1:end-1);
[~,T_max_descent]=findpeaks(-dLFP_dt_minus,'minpeakdistance',100);T_max_descent = T_max_descent/Fs;
%%
T_n{1} = T_neg;
[Cor_n,lag,CI_n]=CrossCorrelation(Ts(chan),T_n,'Fs_ts',2000);
f_target = gca;
T_m{1} = T_max_descent;
[Cor_m,~,CI_m]=CrossCorrelation(Ts(chan),T_m,'Fs_ts',2000);
f_source = get(gca,'Children');
copyobj(f_source,f_target);
F = get(f_target,'Children');
%%
F(1).FaceColor = c2;
F(1).EdgeAlpha=0;
F(2).Color = c2;
F(3).FaceColor = c1;
F(3).EdgeAlpha=0;
F(4).Color = c1;
close;
legend('Negative-peak cross-correlation','95% confidence interval','Maximal-descent cross-correlation','95% confidence interval')
xlim([-0.04 0.04]);

%% 2-3 & 2-4

threshold = 8; % Hz
ch = map(map>0);
Ts_all = cell2mat(Ts(ch)');
Tv = 0:0.01:450;
[R,T]=Ts2Rate(Ts_all,Tv);
R = R/numel(ch);
figure;plot(T,R);
[~,p] = findpeaks(R,'minpeakheight',threshold,'minpeakdistance',10);
%%
T_burst = Tv(p);
% C5 - sub-select ictal discharges C5 20~80, 170-243, 320~419.5
% Nu - sub-select ictal discharges Nu 26~52, 177~213
T_burst1 = union(T_burst(T_burst>20 & T_burst<80),T_burst(T_burst>170 & T_burst<243));
T_burst = union(T_burst,T_burst(T_burst>320 & T_burst<419.5));

disp(['Total: ' num2str(numel(T_burst)) ' episodes are detected']);
% LFP_1 = LFP(chan,:);
Tv_LFP = (1:size(LFP,2))/1000;
[~,LOCS] = ismembertol(T_burst,Tv_LFP,10^-6);
% 
%%
for P = 1:numel(LOCS)

%% Sample, please choose Nu P = 120 % make p_center +8
p_center = LOCS(P)+8;
disp(['The episode happened at ' num2str(T_burst(P)) ' second']);
p_v_all = (p_center - 1000):(p_center + 1000);
p_v = 950:1050;
T_v = (p_v-mean(p_v))/Fs;

figure(2)
NR = [0.1 0.5 1]; % % noise
for j = 1:numel(NR)
A = std(LFP(:,p_v_all),0,2);
n_trial = 1000;
Rec = nan(numel(p_v),n_trial);
T_neg = nan(1,n_trial);
T_des = nan(1,n_trial);
for i = 1:n_trial
    LFP_S = filtfilt(b,a,A(chan)*NR(j)*pinknoise(numel(p_v_all))' + LFP(chan,p_v_all)');
    LFP_S = LFP_S(p_v);   
    [~,p_neg] = min(LFP_S);T_neg(i) = T_v(p_neg);
    [~,p_des] = min(diff(LFP_S));T_des(i) = T_v(p_des);
    Rec(:,i) = LFP_S;
end
if j == 1
subplot(4,1,1);
plot(T_v,Rec,'Color',c3);hold on;
xlim([-0.05 0.05]);
ylabel('micro V')
end
% Histogram
subplot(4,1,j+1);
T_all = [T_neg;T_des]';
hist(T_all,-0.05:0.001:0.05);hold on;
xlim([-0.05 0.05]);
ylim([0,n_trial]);
text(-0.04,max(get(gca,'Ylim'))/2,['A_n_o_i_s_e/A_s_i_g_n_a_l = ' num2str(NR(j))]);
set(gca,'YTickLabel',get(gca,'YTick')/1000);
ylabel('Probability')


legend({'negative peak','maximal descent'},'Location','best')
H = findobj(gca,'type','Patch');
H(1).FaceColor = c2;
H(2).FaceColor = c1;
H(1).EdgeColor = c2;
H(2).EdgeColor = c1;
drawnow;

end
xlabel('Second');


%% All episode
NR_list = 0.01*2.^(0:7); % % noise
M = numel(NR_list);
LFP_S = LFP(ch,p_v_all)';
A = std(LFP_S);
N_test = 200;
N_channel = numel(A);
Record = nan(N_test,N_channel,M,2);
for m = 1:M
    tic;
    disp(['Now calculating NR = ' num2str(NR_list(m))]);
    NR = NR_list(m);
    for i = 1:N_test
        % Noise
        S = bsxfun(@times,NR*A,pinknoise(size(LFP_S,1))') + LFP_S;
        % filter
        S = filtfilt(b,a,S);
        % find the result
        S = S(p_v,:);  
        for j = 1:N_channel
            [~,Record(i,j,m,1)] = min(S(:,j)); % neg-peak
            [~,Record(i,j,m,2)] = min(diff(S(:,j))); % max descend
        end
    end
end
Record = Record / Fs;

STD_info = std(Record,0,1);
STD_info_mean = mean(STD_info,2);
STD{P} = squeeze(STD_info_mean);
figure(3);drawnow
plot(NR_list,STD{P});drawnow;
legend({'Negative peak','Maximal descent'},'Location','best');close;

end
%% plot median and quantile distribution
figure(5);
STD2 = cell2mat(STD);
STD2 =reshape(STD2,numel(NR_list),2,[]);
errorbar(NR_list,median(squeeze(STD2(:,1,:)),2),quantile(squeeze(STD2(:,1,:)),0.25,2),quantile(squeeze(STD2(:,1,:)),0.75,2),'.','MarkerSize',20);hold on
errorbar(NR_list,median(squeeze(STD2(:,2,:)),2),quantile(squeeze(STD2(:,2,:)),0.25,2),quantile(squeeze(STD2(:,2,:)),0.75,2),'.','MarkerSize',20);hold on
set(gca,'XScale','log','YScale','log');
legend({'negative peak','maximal descent'})
xlabel('Noise amplitude / Signal amplitude');
ylabel('millisecond')
% For Nu data, just change it to Nu

%% Figure 3
%3-1 & 2
SpatialAnalysisLFP(LFP, map,'T',[0 100],'Fs',1000,'mode','descent','switch_monitor',1,'f_pass',[1 50],'T_burst',70.67)
%% 3-3 & 3-5
Amp = hilbert(filtfilt(b,a,LFP(ch,:)')');
Amp = abs(Amp);
L(1) = plot(Tv_LFP,mean(Amp));hold on;
L(2) = plot(Tv_LFP,quantile(Amp,0.25));L(2).Color = L(1).Color + (1-L(1).Color)*0.8;
L(3) = plot(Tv_LFP,quantile(Amp,0.75));L(3).Color = L(1).Color + (1-L(1).Color)*0.8;
L(4) = plot(Tv_LFP,mean(Amp));L(4).Color = L(1).Color;
xlabel('Second');ylabel('micro V');title('Patient C5 sample seizure');
legend({'Mean low frequency amplitude','25% percentile','75% percentile'},'Location','best');
xlim([0, 100]);
%% 3-4 & 3-6
LFP_local = LFP(:,320*1000:420*1000);
Rec = SpatialAnalysisLFP(LFP_local, map,'T',[0 100],'Fs',1000,'mode',{'neg-peak','descent'},'switch_monitor',0,'f_pass',[1 50],'T_burst',T_burst(T_burst>320 & T_burst<420)-320);
%% 3-7 & 3-8
Rec_LFP = SpatialAnalysisLFP(LFP, map,'Fs',1000,'mode',{'neg-peak','descent'},'switch_monitor',0,'f_pass',[1 50],'T_burst',T_burst);
Rec_MUA = SpatialAnalysis(Ts, map,'switch_plot',0,'T_burst',T_burst);
[ V_err, D_err, H_err ] = SpatialError( Rec_MUA, Rec_LFP );
%% 4-1,2 
Ts_total = cell2mat(Ts(ch)');
[b,a] = fir1(90,[80 150]/(Fs/2));
LFP_local = filtfilt(b,a,LFP(chan,:));
plot(T_v,abs(hilbert(LFP_local)));hold on;plot(T_v,LFP_local);
T_Raster = [Ts{chan},Ts{chan},Ts{chan}]';T_Raster = T_Raster(:);
Raster = repmat([100;200;nan],[numel(Ts{chan}),1]);
plot(T_Raster,Raster,'k');
%%
Gamma_avg = mean(abs(filtfilt(b,a,LFP')'));
[Output,lag,CI] = CrossCorrelation({Ts_total},Gamma_avg,'Fs',Fs,'Fs_ts',10000,'switch_plot',1);
%% 4-3,4
K1=SpatialAnalysisLFP(LFP,map,'mode','weighted','T',[0 70],'f_pass',[80 150],'T_burst',T_burst(200),'switch_plot',1,'Lossfun','L1')
h=colorbar;ylabel(h,'Power, microV^2');caxis([0 5000]);norm(K1.V)
K2=SpatialAnalysisLFP(LFP,map,'mode','weighted','T',[0 70],'f_pass',[80 150],'T_burst',T_burst(200),'switch_plot',1,'Lossfun','L2')
h=colorbar;ylabel(h,'Power, microV^2');caxis([0 5000]);norm(K2.V)

%% 4-5
Rec_Gamma = SpatialAnalysisLFP(LFP,map,'mode','weighted','f_pass',[80 150],'T_burst',T_burst,'switch_plot',0,'Lossfun','L1');
% Speed > 300 -> set to 0
Rec_Gamma.Significance(sum((Rec_Gamma.V).^2,2)>90000) = 0;
Rec_Gamma2 = SpatialAnalysisLFP(LFP,map,'mode','weighted','f_pass',[80 150],'T_burst',T_burst,'switch_plot',0,'Lossfun','L2');
% Speed > 300 -> set to 0
    % Rotate the result if necessary
    Theta = p.Results.theta;
    V = V * [cos(Theta) sin(Theta);sin(-Theta) cos(Theta)];Rec_Gamma2.Significance(sum((Rec_Gamma2.V).^2,2)>90000) = 0;
%%

V = Rec_Gamma2.V(582+(1:447),:);H = Rec_Gamma2.Significance(582+(1:447));
f_3 = compass(V(H,1) + V(H,2)*sqrt(-1));
T_vector_2 = T_burst(582+(1:447));T_vector_2 = T_vector_2(H);
Color_vector =  (T_vector_2 - min(T_vector_2)) / (max(T_vector_2) - min(T_vector_2));
if length(T_vector_2) > 1              
    for i = 1:numel(T_vector_2)
        f_3(i).Color = [Color_vector(i),0,1-Color_vector(i)];
    end
end
title('Velocity');



%% 4-7, 8
LFP_local = LFP;
Rec_full(1:2) = SpatialAnalysisLFP(LFP_local,map,'channel_index',[],'f_pass',[50],'mode',{'neg-peak','descent'},'T_burst',T_burst,'switch_plot',0,'Lossfun','L2');
Rec_full(3) = SpatialAnalysisLFP(LFP_local,map,'channel_index',[],'f_pass',[80 150],'mode','weighted','T_burst',T_burst,'switch_plot',0,'Lossfun','L2');    
Rec_full(4) = SpatialAnalysisLFP(LFP_local,map,'channel_index',[],'f_pass',[80 150],'mode','weighted','T_burst',T_burst,'switch_plot',0,'Lossfun','L1','Fs',Fs);    
%%
[ V_err, D_err, H_err ] = SpatialError( Rec_MUA, Rec_full(:));    
V_err_full = V_err;
D_err_full = D_err;
H_err_full= H_err;    

%% 5-A-B
T_local = T_burst(884);
p_local = round(T_local*Fs);
p_local = p_local-1000:p_local+1000; %
for i = 1:6
    for j = 1:100 % Remember to turn down shuffle test
        CH = datasample(map_local(map_local>0),(i+2)^2);% use this for 5-1,2
        map_local = map; % use this for 5-1/5-2

        LFP_local = LFP(CH,p_local);
        
        Rec_local_L1(j,1:2) = SpatialAnalysisLFP(LFP_local,map_local,'channel_index',CH,'f_pass',[50],'switch_plot',0,'mode',{'neg-peak','descent'},'T_burst',1,'switch_plot_2',0,'Lossfun','L1');                
        Rec_local_L1(j,3) = SpatialAnalysisLFP(LFP_local,map_local,'channel_index',CH,'f_pass',[80 150],'switch_plot',0,'mode','weighted','T_burst',1,'switch_plot_2',0,'Lossfun','L1');        
        Rec_local_L2(j,1:2) = SpatialAnalysisLFP(LFP_local,map_local,'channel_index',CH,'f_pass',[50],'switch_plot',0,'mode',{'neg-peak','descent'},'T_burst',1,'switch_plot_2',0,'Lossfun','L2');                
        Rec_local_L2(j,3) = SpatialAnalysisLFP(LFP_local,map_local,'channel_index',CH,'f_pass',[80 150],'switch_plot',0,'mode','weighted','T_burst',1,'switch_plot_2',0,'Lossfun','L2');                        
        for m = 1:3
            Vel(i,m,1,j,:) = Rec_local_L1(j,m).V;
            Hel(i,m,1,j,:) = Rec_local_L1(j,m).Significance &  (sum(abs(Rec_local_L1(j,m).V),2)<300) ;                    
            Vel(i,m,2,j,:) = Rec_local_L2(j,m).V;        
            Hel(i,m,2,j,:) = Rec_local_L2(j,m).Significance &  (sum(abs(Rec_local_L2(j,m).V),2)<300);      
        end
    end
end
%% 5A plot
for i = 1:6
    figure;
    for m = 1:3
        L(m)=plot(squeeze(Vel(i,m,2,:,1)),squeeze(Vel(i,m,2,:,2)),'.');  
        [~,~,Lat] = pca(squeeze(Vel(i,m,2,(squeeze(Hel(i,m,2,:,:))),:)));
        % [~,~,Lat] = pca(squeeze(Vel(i,m,2,(squeeze(Hel(i,m,2,:,:))),:)));        
        VAR(i,m) = sum(Lat);
        hold on;
        n_sig(i,m) = sum(logical(squeeze(Hel(i,m,2,:,:))))/100;
    end
    L(4)=plot(squeeze(Vel(i,3,1,:,1)),squeeze(Vel(i,3,1,:,2)),'.');  
    [~,~,Lat] = pca(squeeze(Vel(i,3,1,logical(squeeze(Hel(i,3,1,:,:))),:)));    
    % [~,~,Lat] = pca(squeeze(Vel(i,3,1,:,:)));      
    VAR(i,4) = sum(Lat);    
    n_sig(i,4) = sum(logical(squeeze(Hel(i,3,1,:,:))))/100;    
    legend({['Negative peak'],['Maximal descent'],['L2 weighted'],['L1 weighted']},'location','best');
    % title([num2str(i*10) '% natural background noise']);
    title([num2str((i+2)^2) ' channels used for estimation']);    
    xlabel('cm/sec');ylabel('cm/sec');
    L(1).Color=c1;L(2).Color=c2;L(3).Color=c3;L(4).Color=c4;    
end
%% 5B
figure;
L = plot(VAR(1:6,:)); % For fig 5-1,2
% L = plot(n_sig); % fig 6-1, 2,5,6
set(gca,'XTick',1:6);
set(gca,'XTickLabel',(3:8).^2);
% ylabel('Percentage of successful detection ');
ylabel('Variance of estimator, (cm/sec)^2')
xlabel('Available channels');
legend({['Negative peak'],['Maximal descent'],['L2 weighted'],['L1 weighted']},'location','best');


%% 5C-E
for i = 1:6
    for j = 1:numel(T_burst)
        CH = randsample(map(map>0),(i+2)^2);
        n_pt = round(T_burst(j)*1000);
        LFP_local = LFP(CH,n_pt-1000:n_pt+1000);
        Rec_low_electrode(i,1:2,j) = SpatialAnalysisLFP(LFP_local,map,'channel_index',CH,'f_pass',[1 50],'mode',{'neg-peak','descent'},'T_burst',1,'switch_plot',0,'Lossfun','L2');
        Rec_low_electrode(i,3,j) = SpatialAnalysisLFP(LFP_local,map,'channel_index',CH,'f_pass',[80 150],'mode','weighted','T_burst',1,'switch_plot',0,'Lossfun','L2');    
        Rec_low_electrode(i,4,j) = SpatialAnalysisLFP(LFP_local,map,'channel_index',CH,'f_pass',[80 150],'mode','weighted','T_burst',1,'switch_plot',0,'Lossfun','L1');    
    end
end

%% Fit it in the format
Rec_LE(7,:) = Rec_full; % Add the full electrode case 
for i = 1:7
    if i < 7
    for j = 1:4
        Rec_LE(i,j).T = ([Rec_low_electrode(i,j,:).T]  + T_burst - 1)';
        Rec_LE(i,j).V = reshape([Rec_low_electrode(i,j,:).V],2,[])';
        Rec_LE(i,j).direction = [Rec_low_electrode(i,j,:).direction]';
        Rec_LE(i,j).Significance = [Rec_low_electrode(i,j,:).Significance]';       
        Rec_LE(i,j).method = {Rec_low_electrode(i,j,:).method}';         
        Rec_LE(i,j).pValue = [Rec_low_electrode(i,j,:).pValue]';                 
    end
    end
    [ V_err, D_err, H_err ] = SpatialError( Rec_MUA, Rec_LE(i,:));  
    V_err_lowE{i} = V_err;
    D_err_lowE{i} = D_err;
    H_err_lowE{i} = H_err;    
    close all;    
end
0.65 0.59 0.59 0.62
6.8  13.5  13.8  10.4
%% 5C % of ictal traveling waves as a function of # of electrodes
x_grid = [(3:8).^2,numel(ch)];
L = reshape(cellfun(@sum,{Rec_LE.Significance}),7,[])/numel(T_burst);
L = plot(x_grid,L);
for i = 1:4
    L(i).Color = COLOR(i,:);
    L(i).Marker = 'p';    
    L(i).MarkerFaceColor = COLOR(i,:);    
    L(i).LineStyle = '--';
end
legend({'Neg-peak','Max descent','L2-weighted','L1 weighted'},'Location','Best')
ylabel('% discharges detected as traveling waves');xlabel('number of electrodes available')
ylim([0 1])

%% 5D
for i = 1:7
    Temp = struct2table(Rec_LE(i,:));
    H1 = cell2mat(Temp.Significance');
    H0 = Rec_MUA.Significance(:);
    HH = bsxfun(@and,H0,H1);
    for j = 1:4
        V_ERR(i,j) = median(V_err_lowE{i}(HH(:,j),j));
        V_ERR_U(i,j) = quantile(V_err_lowE{i}(HH(:,j),j),0.75);
        V_ERR_D(i,j) = quantile(V_err_lowE{i}(HH(:,j),j),0.25);
        D_ERR_STD(i,j) = circ_std(D_err_lowE{i}(HH(:,j),j));
        P_ERR(i,j) = circ_rtest(D_err_lowE{i}(HH(:,j),j));
    end
end

x_grid = [(3:8).^2,numel(ch)];
L=errorbar(repmat(x_grid',[1 4]),V_ERR,V_ERR_D,V_ERR_U);
for i = 1:4
    L(i).Color = COLOR(i,:);
    L(i).Marker = 'p';    
    L(i).MarkerFaceColor = COLOR(i,:);    
    L(i).LineStyle = '--';
end
legend({'Neg-peak','Max descent','L2-weighted','L1 weighted'},'Location','Best')
title('Speed error');ylabel('cm/sec');xlabel('# of available electrodes')
%% 5E
L=plot(x_grid,D_ERR_STD);
for i = 1:4
    L(i).Color = COLOR(i,:);
    L(i).Marker = 'p';    
    L(i).MarkerFaceColor = COLOR(i,:);    
    L(i).LineStyle = '--';  
end
legend({'Neg-peak','Max descent','L2-weighted','L1 weighted'},'Location','Best')
title('Direction error');ylabel('radian');xlabel('# of available electrodes')
ylim([0 pi/2]);

%% 6A-B
T_local = T_burst(884);
p_local = round(T_local*Fs);
p_local = p_local-1000:p_local+1000; %
LFP_local = LFP(:,p_local);
for i = 1:11
    for j = 1:100 % Remember to turn down shuffle test
        CH = ch;% use this for 6-1,2
        map_local = map; 

        LFP_local = LFP(CH,p_local);
        for k = 1:size(LFP_local,1)
            NOISE(k,:) = 10 * (i-1) * pinknoise(size(LFP_local,2));
        end
        LFP_local = LFP_local + NOISE;
        
        Rec_local_L1(j,1:2) = SpatialAnalysisLFP(LFP_local,map_local,'channel_index',CH,'f_pass',[50],'switch_plot',0,'mode',{'neg-peak','descent'},'T_burst',1,'switch_plot_2',0,'Lossfun','L1');                
        Rec_local_L1(j,3) = SpatialAnalysisLFP(LFP_local,map_local,'channel_index',CH,'f_pass',[80 150],'switch_plot',0,'mode','weighted','T_burst',1,'switch_plot_2',0,'Lossfun','L1');        
        Rec_local_L2(j,1:2) = SpatialAnalysisLFP(LFP_local,map_local,'channel_index',CH,'f_pass',[50],'switch_plot',0,'mode',{'neg-peak','descent'},'T_burst',1,'switch_plot_2',0,'Lossfun','L2');                
        Rec_local_L2(j,3) = SpatialAnalysisLFP(LFP_local,map_local,'channel_index',CH,'f_pass',[80 150],'switch_plot',0,'mode','weighted','T_burst',1,'switch_plot_2',0,'Lossfun','L2');                        
        for m = 1:3
            Vel(i,m,1,j,:) = Rec_local_L1(j,m).V;
            Hel(i,m,1,j,:) = Rec_local_L1(j,m).Significance &  (sum(abs(Rec_local_L1(j,m).V),2)<300) ;                    
            Vel(i,m,2,j,:) = Rec_local_L2(j,m).V;        
            Hel(i,m,2,j,:) = Rec_local_L2(j,m).Significance &  (sum(abs(Rec_local_L2(j,m).V),2)<300);      
        end
    end
end
%% 6A plot
for i = 1:11
    figure;
    for m = 1:3
        L(m)=plot(squeeze(Vel(i,m,2,:,1)),squeeze(Vel(i,m,2,:,2)),'.');  
        [~,~,Lat] = pca(squeeze(Vel(i,m,2,(squeeze(Hel(i,m,2,:,:))),:)));
        % [~,~,Lat] = pca(squeeze(Vel(i,m,2,(squeeze(Hel(i,m,2,:,:))),:)));        
        VAR(i,m) = sum(Lat);
        hold on;
        n_sig(i,m) = sum(logical(squeeze(Hel(i,m,2,:,:))))/100;
    end
    L(4)=plot(squeeze(Vel(i,3,1,:,1)),squeeze(Vel(i,3,1,:,2)),'.');  
    [~,~,Lat] = pca(squeeze(Vel(i,3,1,logical(squeeze(Hel(i,3,1,:,:))),:)));    
    % [~,~,Lat] = pca(squeeze(Vel(i,3,1,:,:)));      
    VAR(i,4) = sum(Lat);    
    n_sig(i,4) = sum(logical(squeeze(Hel(i,3,1,:,:))))/100;    
    legend({['Negative peak'],['Maximal descent'],['L2 weighted'],['L1 weighted']},'location','best');
    % title([num2str(i*10) '% natural background noise']);
    title([num2str((i-1)*10) ' microV background noise']);    
    xlabel('cm/sec');ylabel('cm/sec');
    L(1).Color=c1;L(2).Color=c2;L(3).Color=c3;L(4).Color=c4;    
    xlim([-100 100]);ylim([0 100])    
end
%% 6B plot
figure;
L = semilogy(10:10:100,VAR(2:end,:)); % For 
% L = plot(n_sig); % fig 6-1, 2,5,6
% ylabel('Percentage of successful detection ');
ylabel('Variance of estimator, (cm/sec)^2')
xlabel('Background noise (microV)');
legend({['Negative peak'],['Maximal descent'],['L2 weighted'],['L1 weighted']},'location','best');

%% 6-C-E background noise simulation
Rec_background(1,:) = Rec_full;
for i = 1:11
    LFP_local = LFP;
    for k = 1:size(LFP_local,1)
        NOISE(k,:) = (i-1) * 10 * pinknoise(size(LFP_local,2));
    end
    LFP_local = LFP_local + NOISE;
    Rec_background(i,1:2) = SpatialAnalysisLFP(LFP_local,map,'f_pass',[50],'mode',{'neg-peak','descent'},'T_burst',T_burst,'switch_plot',0,'Lossfun','L2');
    Rec_background(i,3) = SpatialAnalysisLFP(LFP_local,map,'f_pass',[80 150],'mode','weighted','T_burst',T_burst,'switch_plot',0,'Lossfun','L2');    
    Rec_background(i,4) = SpatialAnalysisLFP(LFP_local,map,'f_pass',[80 150],'mode','weighted','T_burst',T_burst,'switch_plot',0,'Lossfun','L1');             
    [ V_err, D_err, H_err ] = SpatialError( Rec_MUA, Rec_background(i,:));  
    V_err_noise{i} = V_err;
    D_err_noise{i} = D_err;
    H_err_noise{i} = H_err;    
    close all;     
end

%% 6C % of ictal traveling waves as a function of noise
x_grid = 0:10:100;
L = reshape(cellfun(@sum,{Rec_background.Significance}),11,[])/numel(T_burst);
L = plot(x_grid,L);
for i = 1:4
    L(i).Color = COLOR(i,:);
    L(i).Marker = 'p';    
    L(i).MarkerFaceColor = COLOR(i,:);    
    L(i).LineStyle = '--';
end
legend({'Neg-peak','Max descent','L2-weighted','L1 weighted'},'Location','Best')
ylabel('% discharges detected as traveling waves');xlabel('Bbackground noise (microV)')
ylim([0 1])
%% 6D
for i = 1:11
    Temp = struct2table(Rec_background(i,:));
    H1 = cell2mat(Temp.Significance');
    H0 = Rec_MUA.Significance(:);
    HH = bsxfun(@and,H0,H1);
    for j = 1:4
        V_ERR(i,j) = median(V_err_noise{i}(HH(:,j),j));
        V_ERR_U(i,j) = quantile(V_err_noise{i}(HH(:,j),j),0.75);
        V_ERR_D(i,j) = quantile(V_err_noise{i}(HH(:,j),j),0.25);
        D_ERR_STD(i,j) = circ_std(D_err_noise{i}(HH(:,j),j));
        P_ERR(i,j) = circ_rtest(D_err_noise{i}(HH(:,j),j));
    end

end
L=errorbar(repmat(x_grid',[1 4]),V_ERR,V_ERR_D,V_ERR_U);
for i = 1:4
    L(i).Color = COLOR(i,:);
    L(i).Marker = 'p';    
    L(i).MarkerFaceColor = COLOR(i,:); 
    L(i).LineStyle = '--';
end
legend({'Neg-peak','Max descent','L2-weighted','L1 weighted'},'Location','Best')
title('Speed error');ylabel('cm/sec');xlabel('Background noise (microV)')
%% 6E
L=plot(x_grid,D_ERR_STD);
for i = 1:4
    L(i).Color = COLOR(i,:);
    L(i).Marker = 'p';    
    L(i).MarkerFaceColor = COLOR(i,:);    
    L(i).LineStyle = '--';    
end
legend({'Neg-peak','Max descent','L2-weighted','L1 weighted'},'Location','Best')
title('Direction error');ylabel('radian');xlabel('Bbackground noise (microV)')
ylim([0 pi/2]);
%% Figure 1
% C5 unstable/dead electrodes, 
% 1-2, 1-3
load('/mnt/mfs/home/NIMASTER/jl4115/Epilepsy Codes/data/Human data/method paper data/C5.mat')
%%
SpatialAnalysis(Ts,electrodepinout,'T',[350 360],'switch_monitor',0,'T_analysis',[355 356],'min_burst',10)
%%
l = legend('MUA','Speed = 47.57 cm/sec');set(l,'Location','best');
figure(2);
l = legend('MUA','Regression');
set(l,'Location','best');

%% McNemar ... etc.
CC = crosstab(h(:,1),h(:,5))
[h,p,e1,e2] = testcholdout(h(:,1),h(:,2),h(:,5))

%% Make figures pretty
FIG = dir('*.fig');
for i = 7:numel(FIG)
    h = openfig(FIG(i).name);
    h.Position = [1074 369 464 384]
    pubgraph(h,18,4,'w');
    saveas(h,FIG(i).name(1:end-4),'png')
end

%% Figure 7 ECoG
LFP = xltek_eeg.eeg_data;

LFP = bsxfun(@minus,LFP,mean(LFP,2));% Got to de-mean the data
[b,a] = fir1(512,0.5/250,'high');
LFP = filtfilt(b,a,LFP')';
Fs = 500;
T_vector = (1:size(LFP,2))/Fs;
% Plot
plot(T_vector,bsxfun(@plus,LFP(map(1:3,1:3),:),(1:9)'*200));
map = electrodepinout;
interelectrode_distance = 1;
theta = 2*pi*120/360;
Rec_ECoG(1,1:2) = SpatialAnalysisLFP(LFP,map(2:3,1:3),'Fs',500,'f_pass',50,'mode',{'neg-peak','descent'},'switch_plot_2',1,'Lossfun','L2','interelectrode_distance',interelectrode_distance,'dT',0.2,'alpha',0.2,'theta',theta,'T_burst',T_burst_low);
%%
Rec_ECoG(1,3) = SpatialAnalysisLFP(LFP,map(2:3,1:3),'Fs',500,'f_pass',[80 150],'mode','weighted','switch_plot_2',1,'Lossfun','L2','interelectrode_distance',interelectrode_distance,'alpha',0.5,'theta',theta,'T_burst',T_burst_high);
%%
Rec_ECoG(1,4) = SpatialAnalysisLFP(LFP,map(2:3,1:3),'Fs',500,'f_pass',[80 150],'mode','weighted','switch_plot_2',1,'Lossfun','L1','interelectrode_distance',interelectrode_distance,'alpha',0.5,'theta',theta,'T_burst',T_burst_high);
%%
ch = map(1:3,1:3);
T_vector = T_vector(100*500:end);
T_vector = T_vector - min(T_vector);
LFP = LFP(ch,100*500:end);
figure;

[b,a] = fir1(90,50/250);
LFP_l = filtfilt(b,a,LFP')';
[b,a,] = fir1(90,[80 150]/250);
LFP_h = filtfilt(b,a,LFP')';
figure;subplot(3,1,1:2);
L = plot(T_vector,LFP);hold on;for i = 1:9;L(i).Color = [0.9 0.9 0.9];end
plot(T_vector,mean(abs(hilbert(LFP_l))),'b');ylabel('micro V');
subplot(3,1,3);
plot(T_vector,mean(abs(hilbert(LFP_h))),'Color',[0.4660    0.6740    0.1880]);ylabel('micro V');xlabel('second')
temp1 = mean(abs(hilbert(LFP_h)));hold on;
plot(T_vector(round(index)),temp1(round(index)),'o')
%%
subplot(3,1,1:2);hold on;
temp1 = mean(abs(hilbert(LFP_l)));hold on;
index = round((Rec_ECoG(1).T-100)*500);
O = plot(T_vector(round(index)),temp1(round(index)),'o');
%% Adjust 7-0 to 7-4
T_burst_high = T_burst_high(:);
Color_code = [T_burst_high-min(T_burst_high),0*T_burst_high,max(T_burst_high)-T_burst_high]/(max(T_burst_high) - min(T_burst_high));
subplot(3,1,3)
K=allchild(gca);

T_burst_low = T_burst_low(:);
Color_code = [T_burst_low-min(T_burst_low),0*T_burst_low,max(T_burst_low)-T_burst_low]/(max(T_burst_low) - min(T_burst_low));
%% New 7-1 to 7-4
for i = 1:4
    Color_code = [Rec_ECoG(i).T(:)-min(Rec_ECoG(i).T(:)),0*Rec_ECoG(i).T(:),max(Rec_ECoG(i).T(:))-Rec_ECoG(i).T(:)]/(max(Rec_ECoG(i).T(:)) - min(Rec_ECoG(i).T(:)));
    Size = (1 - Rec_ECoG(i).pValue)*50+1;
    figure;
    scatter(Rec_ECoG(i).V(:,1),Rec_ECoG(i).V(:,2),Size,Color_code,'filled')  
    xlim([-150 150]);ylim([-150 150]);grid on;
    xlabel('Vx, cm/sec');ylabel('Vy, cm/sec')
end
%%
hold on;Px=-125*ones(4,1);Py=-80:-20:-140;pV=[0.5,0.1,0.05,0.01];Size = (1 - pV)*50+1;scatter(Px,Py,Size,'filled','k');
% Add legend to size
for i = 1:4
    text(Px(i)+20,Py(i),['p=' num2str(pV(i))]);
end
%% Find corresponding burst T in LFP

L = LFP(ch,:);
L = abs(hilbert(LFP')).^2;
T_vector = (1:size(L,1))/1000-300;
L = mean(L,2);
T_burst = T_burst(end-446:end);
%%
T_burst = T_burst(:);
plot(T_vector,L);
for i = 1:numel(T_burst)
    L_local = L(round(T_burst(i) * 1000) + [-100:100]);
    [L_max(i),ID(i)] = max(L_local);
    ID(i) = ID(i) + round(T_burst(i) * 1000) - 100;
    T_burst_LFP(i) = ID(i)/1000 - 300;
end
Color_code = [T_burst-min(T_burst),0*T_burst,max(T_burst)-T_burst]/(max(T_burst) - min(T_burst));
hold on;
scatter(T_burst_LFP,L(ID),20,Color_code);
xlim([0 150]);
set(gca,'YScale','log');
ylabel('Mean power (microV^2)');
xlabel('Second');
% plot([-300 150],[3 3]*10^4);
plot([20 20 nan 117.5 117.5],[1 10^7 nan 1 10^7]);
legend({'Mean LFP power','Corresponding ictal discharges','Seizure on/offset'},'Location','best')