% Record = SpatialAnalysis(Ts, map, varargin)
%
% Detect burst propagation direction by analyzing spike trains 
%
% Input: Ts - An n_channel by 1 cell having timestamps of each channel
%
%        [Optional]
%        map - map matrix, if left empty - default: electrodepinout.m
%
%        [Optional] 'parameter_name', parameter
%        varargins - Lossfun: (L1/L2), the definition of loss function in
%                             regression
%                    T_burst: Episodes that are already defined, this
%                    will make the program directly skip the episode
%                    selection process
%                    MinPeakHeight: the minimal intensity to define a burst
%                    switch_plot: plot individual regression result
%                    switch_plot_2: plot the whole result
%                    Burst_dT: The length analysis of each burst episode
%                    alpha: Significance level
%                    - Only effective while you don't give T_burst-
%                    T: Time range for evaluating firing rate - 
%                    T_analysis: Time range to detect burst & do analysis
%                    dT: minimal distance of two burst episodes
%                                      
% Output: T - Time of the burst
%         V - velocity
%         direction - direction
%         Significance - statistical test result
%         pValue - pValue
% 
% This function calls SpatialLinearRegression
%                     Ts2Rate
%
% Author: Jyun-you Liou
% Final update: 2016/1/15

function Record = SpatialAnalysis(Ts, map, varargin)
p = inputParser;
addParameter(p,'Lossfun','L2'); % The choice of loss function
addParameter(p,'T',[]); % Range of data to be considered
addParameter(p,'T_analysis',[]); % Range of data to be analyzed, empty = GUI input
addParameter(p,'T_burst',[]); % Pre-selected episodes 
addParameter(p,'switch_plot',0); % Control whether each episode will be monitored realtime or not
addParameter(p,'switch_plot_2',0); % Control plotting the whole result or not
addParameter(p,'dT',0.1); % minimal distance of two burst episodes
addParameter(p,'Burst_dT',0.1); % The range of data analysis for one burst episode
addParameter(p,'MinPeakHeight',[]); % if empty, open a GUI to let the user select
addParameter(p,'alpha',0.05); % Significance level of statistical test
addParameter(p,'channel_index', []); % The data from Ts{i} will be mapped to channel_index(i)
addParameter(p,'interelectrode_distance',0.04); % unit: cm
addParameter(p,'theta',0); % rotate the result in radian
parse(p,varargin{:});

% Set the electrode map 
if nargin < 2
    map = electrodepinout;
elseif isempty(map)
    map = electrodepinout;
end

% Align the channel index to the timestamps data
if p.Results.channel_index
    [Ch,IA,IB] = intersect(map,p.Results.channel_index); 
else
    [Ch,IA,IB] = intersect(map,1:numel(Ts));
end
map(~ismember(map,Ch)) = -1;
n_channel = numel(Ch);

% Extract the physical coordinate of the map
[n1,n2] = size(map);
[p1,p2] = ind2sub([n1,n2],IA);
P = [p1,p2] * p.Results.interelectrode_distance; % unit: cm

% Extract Ts within the selected range of time
Ts = Ts(IB); % Only analyze those channels
Ts = cellfun(@(x) x(:),Ts,'Uniformoutput',0); % Make sure it is suitable for concatenation
Ts = Ts(:);
Ts_all = cell2mat(Ts);
if isempty(p.Results.T)
    T = [min(Ts_all),max(Ts_all)];
else
    T = p.Results.T;
    Ts = cellfun(@(x) x(and(x>T(1),x<T(2))),Ts,'UniformOutput',0);
end
dt = 0.02; % second
T_vector = T(1):dt:T(2);


if isempty(p.Results.T_burst) % Open GUI for user to define burst episodes
    % Convolution to estimate firing rate
    R = Ts2Rate(Ts_all,T_vector)/n_channel; % Convert the unit back to Hz per channel
    % Pick up burst episodes
    figure('units','normalized','outerposition',[0,0.1,1,0.9]);
    hold on;
    plot(T_vector,R);
    % Choose the range of data to detect burst and do analysis
    if isempty(p.Results.T_analysis) 
        xlabel('Second');
        ylabel('Average firing rate, Hz');
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

    % Detect bursts
    Selector = and(T_vector>T_select(1),T_vector<T_select(2));
    dT = p.Results.dT; 
    [~,LOCS] = findpeaks(R .* Selector,'MinPeakDistance',round(dT/dt),'MinPeakHeight',MinPeakHeight);
    plot(T_vector(LOCS),R(LOCS),'ro');

    % Calculate # of burst & define burst data analysis range 
    T_burst = T_vector(LOCS);    
else
    T_burst = p.Results.T_burst;
end
N_burst = numel(T_burst);
Burst_dT = p.Results.Burst_dT;

% Prepare recording variables
V = nan(N_burst,2);
H = false(N_burst,1);
pValue = nan(N_burst,1);

% Analyze burst by burst
for i = 1:N_burst  
    % Retreive the associate spikes   
    Burst_Ts = cellfun(@(x) x( and(x>(T_burst(i)-Burst_dT/2),x<(T_burst(i)+Burst_dT/2))), Ts, ...
                       'UniformOutput',0);
    % Regression
    [V(i,:), H(i), pValue(i)] = SpatialLinearRegression(Burst_Ts,P, ...
                                                        'alpha',p.Results.alpha, ...
                                                        'switch_plot',p.Results.switch_plot, ...
                                                        'Lossfun',p.Results.Lossfun);    
      
end
% Rotate the result if necessary
Theta = p.Results.theta;
V = V * [cos(Theta) sin(Theta);sin(-Theta) cos(Theta)];

% Calculate the velocity in cm/second
Speed = sqrt(sum(V.^2,2));
direction = angle(V(:,1) + V(:,2)*sqrt(-1));

if p.Results.switch_plot_2
    subplot(2,1,1),plot(T_burst,Speed,'.');hold on;
    subplot(2,1,1),plot(T_burst(H),Speed(H),'r*');
    title('Speed');
    ylabel('cm/second');
    try
        ylim([0 max(Speed(H)*1.2)]);
    catch
    end
    subplot(2,1,2),plot(T_burst,direction,'.');hold on;
    subplot(2,1,2),plot(T_burst(H),direction(H),'r*');
    title('Direction')

    % Polar plot
    % In this case we only plot those pass significance test
    figure;
    T_burst_H = T_burst(H);
    Color_vector =  (T_burst_H - min(T_burst_H)) / (max(T_burst_H) - min(T_burst_H));
    f_3 = compass(V(H,1) + V(H,2)*sqrt(-1));
    if length(T_burst_H) > 1              
        for i = 1:length(T_burst_H)
            f_3(i).Color = [Color_vector(i),0.5,1-Color_vector(i)];
        end
    end
    title('Velocity');
end
% Record the result:
Record.T = T_burst; % timing of the burst
Record.V = V; % Velocity, unit: cm/second
Record.direction = direction; % Direction in radian
Record.Significance = H; % Significance
Record.pValue = pValue; % pValue
end