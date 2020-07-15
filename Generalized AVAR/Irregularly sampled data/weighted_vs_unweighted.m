close all
clear all
clc


%% SETUP 
% number of samples
n = 500;
% simulation duration
duration = 1; % [sec]
% number of window lengths
m = 50;
% number of Monte-Carlo simulations
num_monte = 1;
% noise type: 'Gaussain'/'uniform'/'flicker'
noise_type = 'Gaussian';
% amount of noise (for Gaussian is std; for uniform is BW; for flicker is scaler)
noise_gain = 0;
% time stamp sampling method: 'irregular','regular','clustered'
sampling_method = 'irregular';
%%

% create a random walk signal
data_y = random_walk(1000);
data_t = linspace(0,duration,numel(data_y));

avar_regular = zeros(m,num_monte);
avar_irregular_weighted = zeros(m,num_monte);
avar_irregular_unweighted = zeros(m,num_monte);

for k=1:num_monte

    % uniformly sample time stamps
    t_irreg = duration*sort(rand(1,n));

    % regularly sampled time stamps
    t_reg = linspace(0,duration,n);


    % generate noisy measurements
    y_true_reg = get_truth_at(t_reg,data_t, data_y);
    y_true_irreg = get_truth_at(t_irreg,data_t, data_y);
    
    if isequal(noise_type,'Gaussian')
        y_reg = y_true_reg + normrnd(0,noise_gain,[1 n]); % Gaussian white noise
        y_irreg = y_true_irreg + normrnd(0,noise_gain,[1 n]); % Gaussian white noise
    elseif isequal(noise_type,'uniform')
        y_reg = y_true_reg + noise_gain*(rand(1, n)-0.5); % Gaussian white noise
        y_irreg = y_true_irreg + noise_gain*(rand(1, n)-0.5); % Gaussian white noise
    elseif isequal(noise_type,'flicker')
        y_reg = y_true_reg + noise_gain*flicker(n);   % Flicker noise
        y_irreg = y_true_irreg + noise_gain*flicker(n);   % Flicker noise
    else
        disp('Noise type is not defined!')
    end


t_max = max(t_reg);
t_min = min(t_reg);
t_range = t_max - t_min;

% generate a list of 'm' potential window lengths (exponentially sampled)
gamma = 1.1;
tau = t_range/gamma^12;
for i= 1:m-1
    tau = [tau(1)/gamma tau];
end

% calculate Allan variance of the data
avar_regular(:,k) = AVAR2(t_reg,y_reg,tau);
avar_irregular_weighted(:,k) = AVAR2(t_irreg,y_irreg,tau);
avar_irregular_unweighted(:,k) = AVAR2_weigh_off(t_irreg,y_irreg,tau);
end

experiment_name = 'weighted_vs_unweighted';
plot_results_for_this(mean(avar_regular,2),mean(avar_irregular_weighted,2),mean(avar_irregular_unweighted,2),tau,experiment_name);


function plot_results_for_this(avg_avar_regular,avg_avar_irregular_weighted,avg_avar_irregular_unweighted,tau,experiment_name)
%% visualization
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',12);




fig2 = figure('Position', [400 10 500 300]);

        hold on
%         avg_avar = mean(avar,2)';
%         std_avar = std(avar,0,2)';
        plot(tau, avg_avar_regular,'LineWidth',1.5, 'Color',[0.2 0.2 0.2],'LineStyle','-');
        plot(tau, avg_avar_irregular_weighted,'LineWidth',1.5, 'Color',[0.2 0.2 0.2],'LineStyle','--');
        plot(tau, avg_avar_irregular_unweighted,'LineWidth',1.5, 'Color',[0.2 0.2 0.2],'LineStyle',':');
        
        legend('Regularly Sampled','Weighted Averaging','Unweighted Averaging',...
            'Location','northwest')
        %         patch([tau fliplr(tau)], [avg_avar-std_avar fliplr(avg_avar+std_avar)], [.2 .2 .9], 'EdgeColor', 'none', 'FaceAlpha',.1);
        xlabel('Window length $\tau [s]$')
        ylabel('AVAR $\sigma^2_A$')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xlim([tau(1) tau(end)]);

saveas(fig2,['example_' experiment_name '_avar'],'svg')   
end



function y = get_truth_at(t,data_t, data_y)
    
    y = zeros(1,numel(t));
    for i=1:numel(t)
        [~,indx] = min(abs(data_t-t(i)));
        y(i) = data_y(indx);
        % y = 1.0*t + 1*sin(t/1).^1+2;

    end
end




