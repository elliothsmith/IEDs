% Record = SpatialAnalysisLFP(LFP, map, varargin)
%
% Detect burst propagation direction by analyzing LFP data from Utah 
% Array. 
%
% Input: LFP: LFP data can be 1) n_channel by n_t numeric matrix
%                             2) filepath - this .mat needs to contain an
%                                n_channel by n_t numeric matrix inside
%                             3) empty - this will open an GUI.  You can
%                                use the GUI to select the .mat file that
%                                contains the LFP numeric matrix
%
%        map: a numeric matrix with integer entries that point out the 
%             location of channels, if left empty, it will looks for 
%                       electrodepinout.m
%
%        [Optional]
%
%        mode: method selection:
%              neg-peak/amplitude/descent/weighted
%              For running multiple methods, please store the string in a 
%              cell, e.g. {'neg-peak', 'descent'}
%              default: every one
%        Fs: Sampling frequency (Hz) - default: 1000 Hz
%        f_pass: Band pass range (Hz) - 4th order Butterworth filter 
%                default - 1 ~ 50 Hz band pass
%        alpha: Significance level, default: 0.05
%        Lossfun: (L1/L2) - the choice of loss function for regression
%        Burst_dT: The duration of each burst episode
%        channel_index: set LFP(i,:) come from map(map == channel_index(i))
%        T_burst: timing of the bursts, if left empty, it will launch the 
%                 selection GUI
%
%        - Parameters only effective if T_burst is not given -
%        T: Time range of plotting the result
%        T_onset: Marker of the onset of seizure activity
%        T_analysis: Time range to detect burst & do analysis
%        dT: minimal distance of two burst episodes
%        MinPeakHeight: the minimal intensity to define a burst
%
%        - Parameters related to plotting-
%        switch_plot: plot results for every episode
%        switch_plot_2: plot the final whole results
%        
%        - Parameters only for the mode 'weighted'
%        Weight is calculated by W = (A-offset)^power
%        offset: numeric (default: 0)
%                empty (open an GUI to select a range of recording & use 
%                       its mean as offset)
%        power: for calculating weight, default: 2
%
%        - Others -
%        inter_electrode distance: unit: cm, default: 0.04
%
% Output: The spatial correlation of each episode of burst
%
% This function calls: SpatialLinearRegression.m
%
% Author: Jyun-you Liou
% Final update: 2016/1/14

function Record = SpatialAnalysisLFP(LFP, map, varargin)
p = inputParser; 
addParameter(p,'mode',{'neg-peak', 'amplitude', 'descent', 'weighted'}); 
    % option 1: neg-peak (the most negative deflection)
    % option 2: amplitude (timing of max instant amplitude)
    % option 3: descent (the most negative derivative)
    % option 4: weighted (for gamma band data)
addParameter(p,'Fs',1000); % Sampling frequency (Hz)
addParameter(p,'f_pass',[1 50]); % Band pass range
addParameter(p,'alpha',0.05); % Significance level of statistical test
addParameter(p,'Lossfun','L2'); % The choice of loss function
addParameter(p,'channel_index',[]); % The channel index for rows of LFP data
addParameter(p,'T_burst',[]); % for directly given burst timing
addParameter(p,'T',[]); % The range of time you want to plot in LFP data
addParameter(p,'T_analysis',[]); % Range of data to be analyzed, empty = GUI input
addParameter(p,'dT',0.1); % minimal distance of two burst episodes (sec)
addParameter(p,'Burst_dT',0.1); % The range of data analysis for one burst episode
addParameter(p,'MinPeakHeight',[]); % if empty, open a GUI to let the user select
addParameter(p,'switch_plot',0); % Plot real time array map for each burst
addParameter(p,'switch_plot_2',0); % Plot the whole result
addParameter(p,'offset',0);
addParameter(p,'power',2);
addParameter(p,'interelectrode_distance',0.04); % cm
addParameter(p,'theta',0); % Rotate the result in radian
parse(p,varargin{:});

T = p.Results.T; % T: The range of time you want to analyse
Fs = p.Results.Fs; % Sampling frequency

% This allows user to either directly gives the whole data matrix or
% just give the filepath for the data or select data from GUI
if nargin < 1
    % If the user doesn't give LFP data, ask the user where the data is    
    [filename, pathname]=uigetfile({'.mat'});
    data_source = matfile([pathname filename]);
elseif ischar(LFP)
    % For the user who gives filepath
    data_source = matfile(LFP);
elseif and(isnumeric(LFP),ismatrix(LFP))
    % For the user who gives 
    data_source = [];
else
    error('Please give proper data.');    
end

% Get LFP matrix if the input for LFP data is filepath or empty
if ~isempty(data_source)
    % If we need to retrieve the data by the file, screen its variables
    % and find the numerical array, which is supposed to the LFP data
    Vars = whos(data_source); % Retrieve the variable names
    for i = 1:length(Vars) 
        if strcmp(Vars(i).class,'double') || strcmps(Var(i).class,'single')
            % According to whether T is given or not, access the proper
            % range of data
            if any(T)
                Pt = T * Fs;
                Pt(1) = Pt(1) + 1;
                LFP = eval(['data_source.' Vars(i).name(Pt(1):Pt(2))]);
            else
                LFP = eval(['data_source.' Vars(i).name]);
            end
        end
    end    
end

% Find the map & Select electrodes that are going to be analyzed
if nargin < 2
    map = electrodepinout;
elseif isempty(map)
    map = electrodepinout;
end

% Check the channel index for each row of LFP data
if p.Results.channel_index
    [Ch,IA,IB] = intersect(map,p.Results.channel_index); 
else
    [Ch,IA,IB] = intersect(map,1:size(LFP,1));
end
map(~ismember(map,Ch)) = -1;

% Extract the physical coordinate of the map
[n1,n2] = size(map);
[p1,p2] = ind2sub([n1,n2],IA);
P = [p1,p2] * p.Results.interelectrode_distance; % unit: cm

%% Process LFP data 
% Only process channels that 
% 1) both appear in map and LFP matrix data
% 2) within the time range selected by the user
LFP = LFP(IB,:);
[n_channel,n_t] = size(LFP);
Fs = p.Results.Fs;
p_vector = 1:n_t;
T_vector = p_vector / Fs; 
if p.Results.T 
    LFP = LFP(:,and((T_vector>min(p.Results.T)), (T_vector<max(p.Results.T))));
    T_vector = T_vector(:,and((T_vector>min(p.Results.T)), (T_vector<max(p.Results.T))));
end
% Band pass the selected LFP data
[b,a] = fir1(90,p.Results.f_pass/(Fs/2));
LFP = filtfilt(b,a,LFP')';
LFP_Hilbert = hilbert(LFP')';
LFP_amplitude = abs(LFP_Hilbert);

if isempty(p.Results.T_burst)
    %% GUI for user to define bursts
    % The upper panel plots the raw trace from every channel
    figure('units','normalized','outerposition',[0,0.15,1,0.85]);
    subplot(3,1,1:2);
    hold on;
    LFP_display = bsxfun(@plus,max(LFP_amplitude(:))*(1:n_channel)/2,LFP');
    plot(T_vector,LFP_display);
    title('Filtered LFP')
    % The lower panel plots instant LFP amplitude
    subplot(3,1,3);
    LFP_amp_mean = mean(LFP_amplitude);
    f1_2 = plot(T_vector,mean(LFP_amplitude));
    hold on;
    temp = get(f1_2,'Color');
    plot(T_vector,quantile(LFP_amplitude,0.75),'Color',(temp+2)/3,'Linestyle','--');
    plot(T_vector,quantile(LFP_amplitude,0.25),'Color',(temp+2)/3,'Linestyle','--');
    title('Amplitude of LFP, please select analysis range');

    % Choose the range of data to detect burst and do analysis
    if isempty(p.Results.T_analysis) 
        title('Please select analysis range');    
        [T_select,~] = ginput(2);
        T_select = [min(T_select) max(T_select)];
    else
        T_select = p.Results.T_analysis;
    end
    plot([T_select(1) T_select(1)],get(gca,'YLim'),'r');
    plot([T_select(2) T_select(2)],get(gca,'YLim'),'r');

    % Decide what is the minimum intensity to be defined as a burst
    if isempty(p.Results.MinPeakHeight)
        title('Please select the minimal intensity of burst definition');
        [~,MinPeakHeight] = ginput(1);
    else
        MinPeakHeight = p.Results.MinPeakHeight;    
    end

    % Find out the LFP amplitude peak
    Sel = and(T_vector>T_select(1),T_vector<T_select(2));
    [~,LOCS,~,~] = findpeaks(LFP_amp_mean .* Sel ,'MinPeakDistance',round(Fs*p.Results.dT),'MinPeakHeight',MinPeakHeight);
    plot(T_vector(LOCS),LFP_amp_mean(LOCS),'ro');
else
    LOCS = round(p.Results.T_burst*Fs);
end
n_burst = length(LOCS);

% Estimate baseline noise level for offset, only useful for 'weighted' 
if isempty(p.Results.offset)
    title('Please select baseline offset estimation region')
    [T_select,~] = ginput(2);
    plot([T_select(1) T_select(1)],get(gca,'YLim'),'g');
    plot([T_select(2) T_select(2)],get(gca,'YLim'),'g');
    drawnow;
    Sel = and(T_vector>T_select(1),T_vector<T_select(2));
    C = mean(LFP_amplitude(:,Sel),2);
else
    C = p.Results.offset;
end

%% Now start to deal with burst one by one
d_point = ceil(Fs*p.Results.Burst_dT); 

%% Apply each method to perform spatial regression
Record = struct;
method_list = cell(0,0);
if ischar(p.Results.mode)
    method_list{1} = p.Results.mode;
else
    method_list = p.Results.mode;
end

for iter = 1:numel(method_list)
    method = method_list{iter};
    disp(['Current method = ' method]);
    
    % Prepare record variables
    V = nan(n_burst,2);
    H = false(n_burst,1);    
    pValue = nan(n_burst,1);

    for i = 1:n_burst
        display(['Analyzing episode ' num2str(i)]);
        % Create the mask to find local peak
        mask = 0*T_vector; % Reset the mask so that t
        mask(LOCS(i)) = 1;
        mask = logical(conv(mask,ones(1,d_point),'same'));
        T_local = T_vector(mask);

        % Perform regression according to the method
        if strcmpi(method,'weighted')
            %% For weighted estimation method            
            rectify = @(x) x.*(x>0);
            F = @(x) x.^(p.Results.power);
            % 3 steps to calculate the weight for linear regression 
            % (Input to the SpatialLinearRegress function)
            A_local = LFP_amplitude(:,logical(mask));
            % Step 1: Remove things below the threshold
            Weight = rectify(bsxfun(@minus, A_local, C)); %
            % Step 2: Apply the shape of the mask
            Weight = bsxfun(@times, Weight, mask(logical(mask)));
            % Step 3: Non-linear transform
            Weight = F(Weight); % cloud weight  
            % Call the regression machine & do shuffle test there
            
            [V(i,:),H(i),pValue(i)] = SpatialLinearRegression(Weight,P, ...
                                      'switch_plot',p.Results.switch_plot, ...
                                      'Fs',p.Results.Fs, ...
                                      'alpha',p.Results.alpha, ...                                      
                                      'Lossfun',p.Results.Lossfun);
            
        elseif any(strcmpi(method,{'neg-peak','pos-peak','amplitude','descent'}))
            %% For single point estimation methods
            switch method
                case 'neg-peak'
                    [~,Pt_burst] = min(LFP(:,mask),[],2);
                case 'pos-peak'
                    [~,Pt_burst] = max(LFP(:,mask),[],2);
                case 'amplitude'
                    [~,Pt_burst] = max(LFP_amplitude(:,mask),[],2);
                case 'descent'
                    [~,Pt_burst] = min(diff(LFP(:,mask),1,2),[],2);
            end
            Burst_T = T_local(Pt_burst);
            
            % Plot LFP according to the Burst_T sorting results            
            if p.Results.switch_plot 
                figure;                
                [~,T_rank] = sort(Burst_T);
                T_rank(isnan(Burst_T)) = [];
                LFP_obj = plot(T_local,LFP(T_rank,mask));
                local_color_vector = linspace(0,1,n_channel);
                for l = 1:length(Ch)
                    LFP_obj(l).Color = [local_color_vector(l),0,1-local_color_vector(l)];
                end
                waitforbuttonpress;
                close;
            end     
            
            % To indicate it is point process to fit in 
            % SpatialLinearRegression, it needs to be changed to cell form  
            Burst_T = num2cell(Burst_T);
            [V(i,:),H(i),pValue(i)]= SpatialLinearRegression(Burst_T,P, ...
                                     'alpha',p.Results.alpha, ...
                                     'switch_plot',p.Results.switch_plot, ...
                                     'Lossfun',p.Results.Lossfun);
        elseif strcmpi(method,{'xcorr'})
            LFP_local = LFP(:,mask);
            [Beta(:,i),pValue(i)]= SpatialLinearRegressionXcorr(LFP_local,P, ...
                                     'switch_plot',p.Results.switch_plot, ...
                                     'Lossfun',p.Results.Lossfun);
                                 
        end
    end
    % Rotate the result if necessary
    Theta = p.Results.theta;
    V = V * [cos(Theta) sin(Theta);sin(-Theta) cos(Theta)];
    % Calculate the velocity in cm/second
    Speed = sqrt(sum(V.^2,2));
    direction = angle(V(:,1) + V(:,2)*sqrt(-1));
    H(Speed>300) = false;
    
    % Plot the whole result
    if p.Results.switch_plot_2
        T_vector_2 = T_vector(LOCS);
        figure;
        subplot(2,1,1),plot(T_vector_2,Speed,'.');hold on;
        subplot(2,1,1),plot(T_vector_2(H),Speed(H),'r*');
        title('Speed');
        ylabel('cm/second');
        try
            ylim([0 max(norm(V(:,H))*1.2)]);
        catch
        end
        subplot(2,1,2),plot(T_vector_2,direction,'.');hold on;
        subplot(2,1,2),plot(T_vector_2(H),direction(H),'r*');
        title('Direction')

        % Polar plot
        % In this case we only plot those pass significance test
        figure;
        T_vector_2 = T_vector(LOCS(H));
        Color_vector =  (T_vector_2 - min(T_vector(LOCS))) / (max(T_vector(LOCS)) - min(T_vector(LOCS)));
        f_3 = compass(V(H,1) + V(H,2)*sqrt(-1));
        if length(T_vector_2) > 1              
            for i = 1:length(T_vector_2)
                f_3(i).Color = [Color_vector(i),0,1-Color_vector(i)];
            end
        end
        title('Velocity');
    end
    
    % Record the result:
    Record(iter).T = T_vector(LOCS); % timing of the burst
    Record(iter).V = V; % Velocity, unit: cm/second
    Record(iter).direction = direction; % Direction in radian
    Record(iter).Significance = H; % Significance
    Record(iter).method = method; % The method to detect discharge at each electrode
    Record(iter).pValue = pValue; % pValue
end

end