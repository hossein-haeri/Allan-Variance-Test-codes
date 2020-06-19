close all
clear all
clc

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',16);


%% SETUP 
% number of samples
n = 200;
% simulation duration
duration = 50; % [sec]
% space length
space_length = 1000; % [sec]
% number of window lengths (time)
m_t = 12;
% number of window lengths (space)
m_x = 1;
% number if Monte-Carlo simulations
num_monte = 1;
% noise type: 'Gaussain'/'uniform'/'flicker'
noise_type = 'Gaussian';
% amount of noise (for Gaussian is std; for uniform is BW; for flicker is scaler)
noise_gain = 1.5;
% time stamp sampling method: 'irregular','regular','clustered'
sampling_method = 'irregular';



%%
avar = zeros(m_t,m_x,num_monte);
mse = zeros(m_t,m_x,num_monte);


for k=1:num_monte
    %% time stamp sampling (regular,irregular uniform, irregular Gaussian)
    if isequal(sampling_method,'irregular')
        % uniformly sample time stamps
%         t = duration*sort(rand(1,n));
%         x = space_length*sort(rand(1,n));
        t = duration*rand(1,n);
        x = space_length*rand(1,n);
    elseif isequal(sampling_method,'regular')
        % regularly sampled time stamps
        t = linspace(0,duration,n);
        x = linspace(0,space_length,n);
    elseif isequal(sampling_method,'clustered')
        % time stamps sampled from a Gaussian distribution
        t = duration*sort(normrnd(0,.2,[1,n]));
        t = t + abs(min(t));
        x = space_length*sort(normrnd(0,.2,[1,n]));
    end
    
    t_max = max(t);
    t_min = min(t);
    t_range = t_max - t_min;
    
    x_max = max(x);
    x_min = min(x);
    x_range = x_max - x_min;
    

    % generate noisy measurements
    y_true = get_truth_at_1D(t,x);
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
    gamma_t = 1.4;
    tau_t = t_range/4;
    for i= 1:m_t-1
            tau_t = [tau_t(1)/gamma_t tau_t];
    end
    
    gamma_x = 1.4;
    tau_x = x_range/4;
    for i= 1:m_x-1
            tau_x = [tau_x(1)/gamma_x tau_x];
    end
    
    


    % calculate Allan variance of the data
    avar(:,:,k) = AVAR3(t,x,y,tau_t,tau_x);
%     avar = AVAR3(t,x,y,tau_t,tau_x);
    % calculate Mean Squared Errror of the moving average estimation
%     mse(:,:,k) = MSE(t,x,y,tau_t,tau_x);

end

figure(1)
surf(avar)

figure(2)
scatter3(t,x,y,60,'.','DisplayName','Records','MarkerEdgeAlpha',0.7)
% %% visualization
%%
% figure(1)
%     hold on
%     T = (t_min:t_range/100:t_max);
%     X = (x_min:x_range/100:x_max);
%     Y = get_truth_at(t,x);
% %     plot(T,get_truth_at(T,X),'k-','LineWidth',2,'DisplayName','Actual value')
%     scatter3(t,x,Y,60,'.','DisplayName','Records','MarkerEdgeAlpha',0.7)
% %     [~,idx] = min(mean(avar,2));
% %     xline(t_range/2);
% %     xline(t_range/2+tau_t(idx));
%     xlabel('$t$')
%     ylabel('$\tilde{\theta}$')
%     legend('Actual value', 'Records')

%%
% figure(2)
%     subplot(2,1,1)
%         hold on
%         avg_avar = mean(avar,2)';
%         std_avar = std(avar,0,2)';
%         ax1 = plot(tau_t, avg_avar,'LineWidth',2);
%         patch([tau_t fliplr(tau_t)], [avg_avar-std_avar fliplr(avg_avar+std_avar)], [.2 .2 .9], 'EdgeColor', 'none', 'FaceAlpha',.2);
%         xlabel('Window length $\tau [s]$')
%         ylabel('AVAR $\sigma^2_\theta$')
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         grid on
%         xl = xlim;
%     subplot(2,1,2)
%         hold on
%         avg_mse = mean(mse,2)';
%         std_mse = std(mse,0,2)';
%         ax2 = plot(tau_t, avg_mse,'LineWidth',2);
%         patch([tau_t fliplr(tau_t)], [avg_mse-std_mse fliplr(avg_mse+std_mse)], [.2 .2 .9], 'EdgeColor', 'none', 'FaceAlpha',.2);
%         xlabel('Window length $\tau [s]$')
%         ylabel('Estimation MSE')
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         grid on
%         xlim(xl);


function L = MSE(t,x,y,tau_t_list,tau_x_list)
L = [];
t_max = max(t);
t_min = min(t);
y_t = get_truth_at((t_min:0.01:t_max))';
initial_guess = 0;
% for each window length tau
for i= 1:numel(tau_x_list)
    for j= 1:numel(tau_t_list)
        y_hat_hist = [];
        tau_t = tau_t_list(j);
        tau_x = tau_x_list(i);
        % start the simulate
        for space=x_min:x_range/1000:x_max
            for time=t_min:t_range/1000:t_max
                % collect data points that fall into the window 
                valid_hist = y(t<=time & t>=time-tau_t, x<=space+tau_x/2 & x>=space-tau_x/2);
                % if window is not empty
                if ~isempty(valid_hist)
                    % estimate the current value as the average over the window
                    y_hat = mean(valid_hist,'all');
                else
                    % for empty window cases
                    y_hat = initial_guess;
                end
                y_hat_hist = [y_hat_hist; y_hat];
            end
        end
        % compute the MSE of the moving average estimate using window length tau
        l = sum((y_t-y_hat_hist).^2./numel(y_t));
        L = [L; l];
    end
end
end



function y = get_truth_at(t,x)
    % change the function for different actual values
    y = zeros(numel(t),numel(x));
    for i_x=1:numel(x)
        for i_t=1:numel(t)
            y(i_t,i_x) = 0.2*t(i_t) + 2*sin(t(i_t)/2) + x(i_x)./50;
        end
    end
end

function y = get_truth_at_1D(t,x)
    % change the function for different actual values
            y = 0.2*t + 2*sin(t/2) + 0*x./50;

end

