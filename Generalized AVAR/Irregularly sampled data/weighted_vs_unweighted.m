close all
clear all
clc


%% SETUP 
% number of samples
n = 1000;
% number of samples of random walk signal
n_RW = 5000;
% simulation duration
duration = 1; % [sec]
% number of window lengths
m = 20;
% number of Monte-Carlo simulations
num_monte = 100;
% noise type: 'Gaussain'/'uniform'/'flicker'
noise_type = 'Gaussian';
% amount of noise (for Gaussian is std; for uniform is BW; for flicker is scaler)
noise_gain_reg = 0;
noise_gain = 0;
% time stamp sampling method: 'irregular','regular','clustered'
sampling_method = 'irregular';
%%

avar_regular = zeros(m,num_monte);
avar_irregular_weighted = zeros(m,num_monte);
avar_irregular_unweighted = zeros(m,num_monte);


for k=1:num_monte
    
    % create a random walk signal
    data_y = random_walk(n_RW);
    data_t = linspace(0,duration,numel(data_y));
    



    % uniformly sample time stamps
    t_irreg = duration*sort(rand(1,n));

    % regularly sampled time stamps
    t_reg = linspace(0,duration,n_RW);


    % generate noisy measurements
    y_true_reg = get_truth_at(t_reg,data_t, data_y);
    y_true_irreg = get_truth_at(t_irreg,data_t, data_y);
    
    if isequal(noise_type,'Gaussian')
        y_reg = y_true_reg + normrnd(0,noise_gain_reg,[1 n_RW]); % Gaussian white noise
        y_irreg = y_true_irreg + normrnd(0,noise_gain,[1 n]); % Gaussian white noise
    elseif isequal(noise_type,'uniform')
        y_reg = y_true_reg + noise_gain_reg*(rand(1, n_RW)-0.5); % Gaussian white noise
        y_irreg = y_true_irreg + noise_gain*(rand(1, n)-0.5); % Gaussian white noise
    elseif isequal(noise_type,'flicker')
        y_reg = y_true_reg + noise_gain_reg*flicker(n_RW);   % Flicker noise
        y_irreg = y_true_irreg + noise_gain*flicker(n);   % Flicker noise
    else
        disp('Noise type is not defined!')
    end


t_max = max(t_reg);
t_min = min(t_reg);
t_range = t_max - t_min;

% generate a list of 'm' potential window lengths (exponentially sampled)
gamma = 1.15;
tau = t_range/gamma^5;
for i= 1:m-1
    tau = [tau(1)/gamma tau];
end

% calculate Allan variance of the data
% avar_regular(:,k) = AVAR(t_reg,y_reg,floor(tau.*n_RW));
avar_regular(:,k) = AVAR(y_reg,floor(tau.*n_RW));
avar_irregular_weighted(:,k) = AVAR2(t_irreg,y_irreg,tau);
avar_irregular_unweighted(:,k) = AVAR2_weigh_off(t_irreg,y_irreg,tau);
end

experiment_name = 'weighted_vs_unweighted';
plot_results_for_this(avar_regular,avar_irregular_weighted,avar_irregular_unweighted,tau,experiment_name);

fig1 = figure('Position', [400 510 500 300]);
w = abs(avar_irregular_weighted-avar_regular);
u = abs(avar_irregular_unweighted-avar_regular);
m_w = mean(w,2);
m_u = mean(u,2);
hold on
plot(tau,m_u,'DisplayName','simple averaging','Marker','s')
plot(tau,m_w,'DisplayName','weighted averaging','Marker','^')
xlabel('Window length $\tau [s]$')
ylabel('AVAR error')
set(gca,'xscale','log')
legend('show','Location','northwest')
grid on
saveas(fig1,'weighted_vs_unweighted_error','svg')


function plot_results_for_this(avg_avar_regular,avg_avar_irregular_weighted,avg_avar_irregular_unweighted,tau,experiment_name)
%% visualization
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',14);




fig2 = figure('Position', [400 10 500 300]);

        hold on
%         avg_avar = mean(avar,2)';
%         std_avar = std(avar,0,2)';

        plot(tau, mean(avg_avar_regular,2),'LineWidth',0.75, 'Color',[0.2 0.2 0.2],'LineStyle','-','DisplayName','regularly sampled');
        plot(tau, mean(avg_avar_irregular_unweighted,2),'LineWidth',0.75, 'Color',[0.2 0.2 0.2],'LineStyle','-.','DisplayName','simple averaging');
        plot(tau, mean(avg_avar_irregular_weighted,2),'LineWidth',0.75, 'Color',[0.2 0.2 0.2],'LineStyle','--','DisplayName','weighted averaging');
%         plot(tau, avg_avar_regular,'LineWidth',1.5, 'Color',[0.2 0.2 0.2],'LineStyle','-');
%         plot(tau, avg_avar_irregular_weighted,'LineWidth',1.5, 'Color',[0.2 0.2 0.2],'LineStyle','--');
%         plot(tau, avg_avar_irregular_unweighted,'LineWidth',1.5, 'Color',[0.2 0.2 0.2],'LineStyle',':');
   
        %         patch([tau fliplr(tau)], [avg_avar-std_avar fliplr(avg_avar+std_avar)], [.2 .2 .9], 'EdgeColor', 'none', 'FaceAlpha',.1);
        xlabel('Window length $\tau [s]$')
        ylabel('AVAR $\sigma^2_A$')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xlim([tau(1) tau(end)]);
        legend('show','Location','northwest')
saveas(fig2,'weighted_vs_unweighted','svg')   
end



function y = get_truth_at(t,data_t, data_y)
    
    y = zeros(1,numel(t));
    for i=1:numel(t)
        [~,indx] = min(abs(data_t-t(i)));
        y(i) = data_y(indx);
        % y = 1.0*t + 1*sin(t/1).^1+2;

    end
end




