close all
clear all
clc

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',16);

load('jerath_signal')

%% SETUP 
% number of samples
n = 1000;
% number of window lengths
m = 30;
% number if Monte-Carlo simulations
num_monte = 10;


avar = zeros(m,num_monte);
mse = zeros(m,num_monte);

for k=1:num_monte

% uniformly sample time stamps
t = 50*sort(rand(1,n));

% % regularly sampled time stamps
% t = linspace(0,10,n);

% % time stamps sampled from a Gaussian distribution
% t = 10*sort(normrnd(0,.2,[1,n]));
% t = t + abs(min(t));

%%



% generate noisy measurements
y_true = get_truth_at(t,tout, simout);
y = y_true + normrnd(0,.001,[1 n]);

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
avar(:,k) = AVAR2(t,y,tau);

% calculate Mean Squared Errror of the moving average estimation
mse(:,k) = MSE(t,y,tau,tout, simout);



end


%% visualization

figure(1)
    hold on
    plot(tout,simout,'k-','LineWidth',1,'DisplayName','Actual value')
    scatter(t,y,60,'.','DisplayName','Records','MarkerEdgeAlpha',0.7)
    [~,idx] = min(mean(avar,2));
    xline(t_range/2);
    xline(t_range/2+tau(idx));
    legend('Actual value', 'Records')
    xlabel('$t$')
    ylabel('$\tilde{\theta}$')
    
figure(2)
    subplot(2,1,1)
        avg_avar = mean(avar,2)';
        std_avar = std(avar,0,2)';
        ax1 = plot(tau, avg_avar,'LineWidth',2);
        patch([tau fliplr(tau)], [avg_avar-std_avar fliplr(avg_avar+std_avar)], [.2 .2 .9], 'EdgeColor', 'none', 'FaceAlpha',.2);
        xlabel('Window length $\tau [s]$')
        ylabel('AVAR $\sigma^2_\theta$')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xl = xlim;
    subplot(2,1,2)
        avg_mse = mean(mse,2)';
        std_mse = std(mse,0,2)';
        ax2 = plot(tau, avg_mse,'LineWidth',2);
        patch([tau fliplr(tau)], [avg_mse-std_mse fliplr(avg_mse+std_mse)], [.2 .2 .9], 'EdgeColor', 'none', 'FaceAlpha',.2);
        xlabel('Window length $\tau [s]$')
        ylabel('Estimation MSE')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        xlim(xl);

        
function L = MSE(t,y,tau_list,tout, simout)
L = [];
t_max = max(t);
t_min = min(t);
y_t = get_truth_at((t_min:0.01:t_max),tout, simout)';
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





function y = get_truth_at(t,tout, simout)
    
    y = zeros(1,numel(t));
    for i=1:numel(t)
        [~,indx] = min(abs(tout-t(i)));
        y(i) = simout(indx);
        % y = 1.0*t + 1*sin(t/1).^1+2;

    end
end
