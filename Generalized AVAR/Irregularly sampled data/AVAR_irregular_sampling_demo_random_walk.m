close all
clear all
clc

rng(115);
%% SETUP 
% number of samples
n = 500;
% simulation duration
duration = 1; % [sec]
% number of window lengths
m = 100;
% number of Monte-Carlo simulations
num_monte = 500;
% noise type: 'Gaussain'/'uniform'/'flicker'
noise_type = 'Gaussian';
% amount of noise (for Gaussian is std; for uniform is BW; for flicker is scaler)
noise_gain = 5;
% time stamp sampling method: 'irregular','regular','clustered'
sampling_method = 'irregular';
%%

% create a random walk signal
data_y = random_walk(5000);
data_t = linspace(0,duration,numel(data_y));

avar = zeros(m,num_monte);
mse = zeros(m,num_monte);

for k=1:num_monte
    %% time stamp sampling (regular,irregular uniform, irregular Gaussian)
    if isequal(sampling_method,'irregular')
        % uniformly sample time stamps
        t = duration*sort(rand(1,n));
    elseif isequal(sampling_method,'regular')
        % regularly sampled time stamps
        t = linspace(0,duration,n);
    elseif isequal(sampling_method,'clustered')
        % time stamps sampled from a Gaussian distribution
        t = duration*sort(normrnd(0,.2,[1,n]));
        t = t + abs(min(t));
    end

    % generate noisy measurements
    y_true = get_truth_at(t,data_t, data_y);
    if isequal(noise_type,'Gaussian')
        y = y_true + normrnd(0,noise_gain,[1 n]); % Gaussian white noise
    elseif isequal(noise_type,'uniform')
        y = y_true + noise_gain*(rand(1, n)-0.5); % Gaussian white noise
    elseif isequal(noise_type,'flicker')
        y = y_true + noise_gain*flicker(n);   % Flicker noise
    else
        disp('Noise type is not defined!')
    end


t_max = max(t);
t_min = min(t);
t_range = t_max - t_min;

% generate a list of 'm' potential window lengths (exponentially sampled)
gamma = 1.05;
tau = t_range/3;
for i= 1:m-1
    tau = [tau(1)/gamma tau];
end

% calculate Allan variance of the data
avar(:,k) = AVAR2(t,y,tau);

% calculate Mean Squared Errror of the moving average estimation
mse(:,k) = MSE(t,y,tau,data_t, data_y);



end
experiment_name = 'random_walk';
analyze_results(t,y,data_t,data_y,avar,mse,tau,experiment_name);
plot_results(t,y,data_t,data_y,avar,mse,tau,experiment_name);


function L = MSE(t,y,tau_list,tout, simout)
L = [];
t_max = max(t);
t_min = min(t);
y_t = get_truth_at((t_min:(t_max-t_min)/100:t_max),tout, simout)';
initial_guess = 0;
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
    l = sum((y_t-y_hat_hist).^2./numel(y_t));
    L = [L; l];
end
end


function y = get_truth_at(t,data_t, data_y)
    
    y = zeros(1,numel(t));
    for i=1:numel(t)
        [~,indx] = min(abs(data_t-t(i)));
        y(i) = data_y(indx);
        % y = 1.0*t + 1*sin(t/1).^1+2;

    end
end
