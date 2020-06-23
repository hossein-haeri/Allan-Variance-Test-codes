close all
clear all
clc

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',16);


%% SETUP 
% number of samples
n = 500;
% simulation time duration
duration = 50; % [sec]
% number of window lengths
m = 30;
% number if Monte-Carlo simulations
num_monte = 1;
% noise type: 'Gaussain'/'uniform'/'flicker'
noise_type = 'uniform';
% amount of noise (for Gaussian is std; for uniform is BW; for flicker is scaler)
noise_gain = 0.5;
% time stamp sampling method: 'irregular','regular','clustered'
sampling_method = 'irregular';



%%
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
    
    t_max = max(t);
    t_min = min(t);
    t_range = t_max - t_min;
    

    % generate noisy measurements
    y_true = get_truth_at(t);
    if isequal(noise_type,'Gaussian')
        y = y_true + normrnd(0,noise_gain,[1 n]); % Gaussian white noise
    elseif isequal(noise_type,'uniform')
        y = y_true + noise_gain*(rand(1, n)-0.5); % Gaussian white noise
    elseif isequal(noise_type,'flicker')
        y = y_true + noise_gain*flicker(n);   % Flicker noise
    else
        disp('Noise type is not defined!')
    end


    % generate a list of 'm' potential window lengths (exponentially sampled)
    gamma = 1.2;
    tau = t_range/4;
    for i= 1:m-1
        tau = [tau(1)/gamma tau];
    end
    T = t(2)-t(1);
    m_list = round(tau./T);
    % calculate Allan variance of the data
    avar(:,k) = AVAR2(tau,y,m_list);

    % calculate Mean Squared Errror of the moving average estimation
    mse(:,k) = MSE(t,y,tau);
end

%% visualization
figure(1)
    hold on
    T = (t_min:t_range/100:t_max);
    plot(T,get_truth_at(T),'k-','LineWidth',2,'DisplayName','Actual value')
    scatter(t,y,60,'.','DisplayName','Records','MarkerEdgeAlpha',0.7)
    [~,idx] = min(mean(avar,2));
    xline(t_range/2);
    xline(t_range/2+tau(idx));
    xlabel('$t$')
    ylabel('$\tilde{\theta}$')
    legend('Actual value', 'Records')
figure(2)
    subplot(2,1,1)
        hold on
        avg_avar = mean(avar,2)';
        std_avar = std(avar,0,2)';
        plot(m_list, avg_avar,'LineWidth',2);
%         patch([m_list fliplr(m_list)], [avg_avar-std_avar fliplr(avg_avar+std_avar)], [.2 .2 .9], 'EdgeColor', 'none', 'FaceAlpha',.2);
        xlabel('Window length $\tau [s]$')
        ylabel('AVAR $\sigma^2_\theta$')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xlim([tau(1) tau(end)]);
    subplot(2,1,2)
        hold on
        avg_mse = mean(mse,2)';
        std_mse = std(mse,0,2)';
        plot(tau, avg_mse,'LineWidth',2);
%         patch([tau fliplr(tau)], [avg_mse-std_mse fliplr(avg_mse+std_mse)], [.2 .2 .9], 'EdgeColor', 'none', 'FaceAlpha',.2);
        xlabel('Window length $\tau [s]$')
        ylabel('Estimation MSE')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xlim([tau(1) tau(end)]);


function L = MSE(t,y,tau_list)
L = [];
t_max = max(t);
t_min = min(t);
y_t = get_truth_at((t_min:0.01:t_max))';
initial_guess = 0;
% for each window length tau
for j= 1:numel(tau_list)
    y_hat_hist = [];
    tau = tau_list(j);
    % start the simulate
    for time=t_min:0.01:t_max
        % collect data points that fall into the window 
        valid_hist = y(t<=time & t>time-tau);
        % if window is not empty
        if ~isempty(valid_hist)
            % estimate the current value as the average over the window
            y_hat = mean(valid_hist);
        else
            % for empty window cases
            y_hat = initial_guess;
        end
        y_hat_hist = [y_hat_hist; y_hat];
    end
    % compute the MSE of the moving average estimate using window length tau
    l = sum((y_t-y_hat_hist).^2./numel(y_t));
    L = [L; l];
end
end



function y = get_truth_at(t)
    % change the function for different actual values
    y = 0.2*t + 1*sin(t/2);
    
end



