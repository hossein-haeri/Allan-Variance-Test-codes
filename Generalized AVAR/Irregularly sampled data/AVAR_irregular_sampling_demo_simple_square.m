close all
clear all
clc



%%%%% SETUP %%%%
% number of samples
n = 1000;

% simulation duration
duration = 1000; % [sec]

% number of window lengths
m = 100;

% number of Monte-Carlo simulations
num_monte = 100;

% noise type: 'Gaussain'/'uniform'/'flicker'
noise_type = 'Gaussian';

% amount of noise
noise_gain = 0.5;  %(for Gaussian is std; for uniform is BW; for flicker is scaler)

% tau sampling setting
gamma = 1.05;
P = 10;

% time stamp sampling method: 'irregular','regular','clustered'
sampling_method = 'irregular';



avar = zeros(m,num_monte);
mse = zeros(m,num_monte);

for k=1:num_monte
   % time stamp sampling (regular,irregular uniform, irregular Gaussian)
    if isequal(sampling_method,'irregular')
        % uniformly sample time stamps
        t = duration*sort(rand(1,n));
    elseif isequal(sampling_method,'regular')
        % regularly sampled time stamps
        t = linspace(0,duration,n);
    elseif isequal(sampling_method,'clustered')
        % time stamps sampled from a Gaussian d5istribution
        t = duration*sort(normrnd(0,.2,[1,n]));
        t = t + abs(min(t));
    end
    
    t_max = max(t);
    t_min = min(t);
    t_range = t_max - t_min;
    tau = t_range/P;
    
    
    % gen 2erate noisy measurements
    y_true = get_truth_at(t);
    if isequal(noise_type,'Gaussian')
        y = y_true + normrnd(0,noise_gain,[1 n]); % Gaussian white noise
    elseif isequal(noise_2type,'uniform')
        y = y_true + noise_gain*(rand(1, n)-0.5); % Gaussian white noise
    elseif isequal(noise_type,'flicker')
        y = y_true + noise_gain*flicker(n);   % Flicker noise
    else
        disp('Noise type is not defined!')
    end

    % generate a list of 'm' potential window lengths (exponentially sampled)

    for i= 1:m-1
        tau = [tau(1)/gamma tau];
    end

    % calculate Allan variance of the data
    avar(:,k) = AVAR2(t,y,tau);

    % calculate Mean Squared Errror of the moving average estimation
    mse(:,k) = MSE(t,y,tau);
    
end

T = (t_min:t_range/200:t_max);
Y = get_truth_at(T);


experiment_name = 'square';

plot_results(t,y,T,Y,avar,mse,tau,experiment_name);
analyze_results(t,y,T,Y,avar,mse,tau,experiment_name);


function L = MSE(t,y,tau_list)
L = [];
t_max = max(t);
t_min = min(t);
y_t = get_truth_at((t_min:(t_max-t_min)/100:t_max))';
y_hat = 0;
% for each window length tau
for j= 1:numel(tau_list)
    y_hat_hist = [];
    tau = tau_list(j);
    % start the simulate
    for time=t_min:(t_max-t_min)/100:t_max
        % collect data points that fall into the window 
        valid_hist = y(t<=time & t>time-tau);
        % if window is not empty
        if ~isempty(valid_hist)
            % estimate the current value as the average over the window
            y_hat = mean(valid_hist);
        end
        y_hat_hist = [y_hat_hist; y_hat];
    end
    % compute the MSE of the moving average estimate using window length tau
    l = sum((y_t-y_hat_hist).^2)./numel(y_t);
    L = [L; l];
end
end



function y = get_truth_at(t)
    % change the function for different actual values
%     y = 1*t;
    y = square(0.02*t);
%     y = sin(2*t)+0;

end

