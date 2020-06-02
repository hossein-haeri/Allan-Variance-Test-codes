close all
clear all
clc

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',16);
% data = readtable('vehicles_data.csv','ReadVariableNames',true);


%% SETUP 
% number of samples
n = 500;
% number of window lengths
m = 30;

% % uniformly sample time stamps
% t = 10*sort(rand(1,n));

% % regularly sampled time stamps
% t = linspace(0,10,n);

% % time stamps sampled from a Gaussian distribution
t = 10*sort(normrnd(0,.2,[1,n]));
t = t + abs(min(t));

%%



% generate noisy measurements
y_true = get_truth_at(t);
y = y_true + normrnd(0,.5,[1 n]);

t_max = max(t);
t_min = min(t);
t_range = t_max - t_min;

% generate a list of 'm' potential window lengths (exponentially sampled)
gamma = 1.2;
tau = t_range/gamma^10;
for i= 1:m-1
    tau = [tau(1)/gamma tau];
end

% calculate Allan variance of the data
avar = AVAR2(t,y,tau);

% calculate Mean Squared Errror of the moving average estimation
mse = MSE(t,y,tau);






%% visualization

figure(1)
    scatter(t,y,'.')
    [~,idx] = min(avar);
    xline(t_range/2);
    xline(t_range/2+tau(idx));
    xlabel('$t$')
    ylabel('$\tilde{\theta}$')
    
figure(2)
    subplot(2,1,1)
        ax1 = plot(tau, avar,'LineWidth',2);
        xlabel('Window length $\tau [s]$')
        ylabel('AVAR $\sigma^2_\theta$')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xl = xlim;
    subplot(2,1,2)
        ax2 = plot(tau, mse,'LineWidth',2);
        xlabel('Window length $\tau [s]$')
        ylabel('Estimation MSE')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xlim(xl);



