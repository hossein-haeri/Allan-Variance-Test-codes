clear all
clc
close all

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',16);

c = lines(4);
simul_timesteps = 50000; % number of the simulation timesteps
num_measurements = 5000; % number of the measurements
SNR = 20; % Signal to Noise power ratio
num_tau = 50; % number of the timescales in AVAR analysis
horizon = .1; % DAVAR horizon [fraction of the entire signal]


% sample timestamps
% timestamps = multi_modal_timestamp_generator(simul_timesteps,num_measurements,5);
timestamps = sort(rand(1,num_measurements));
t_r = linspace(0,1,simul_timesteps);

% build up a reference signal
% [signal, signal_ref] = random_polynomial_generator(t_r, timestamps, 10);
[signal, signal_ref] = random_sine_generator(t_r, timestamps, 50*rand(1,10));


signal = signal - mean(signal);
signal_ref = signal_ref - mean(signal_ref);

% build up a noise signal
noise = normrnd(0,1,[1, num_measurements]);

% amp = 5*sin(20/num_measurements*(1:num_measurements)) + 2;
% noise = noise .* amp;

% regulate signal/noise powers according to the given SNR
% and build the measurement signal as y
signal_power = rms(signal)^2;
noise_power = rms(noise)^2;
snr = signal_power/noise_power;
y = sqrt(SNR/snr)*signal + noise;
signal_ref = signal_ref * sqrt(SNR/snr);

% make an exponentially increasing list of window lengths
tau_max = horizon/2;
tau_min = 2/simul_timesteps;
gamma = (tau_max/tau_min)^(1/num_tau);
tau = tau_min;
for i= 1:num_tau-1
    tau = [tau tau(end)*gamma];
end

h = floor(horizon*simul_timesteps);

y_max = max(y);
y_min = min(y);

subplot(2,1,1)


AVAR = [];
MSE = [];
opt_est = [];
avar_est = [];
tick = 100;
% figure('units','normalized','outerposition',[0 0 1 1])
x_hat = zeros(1,numel(tau));
t_est = [];
ll = 0;
for k= h+1:tick:simul_timesteps
    ll = ll + tick;
    t = k/simul_timesteps;
    recent = timestamps<=t & timestamps > t-horizon;
    recent_r = t_r<=t & t_r > t-horizon;
    
    avar = AVAR2(timestamps(recent),y(recent),tau);
    
    [~,avar_min] = min(avar);
    AVAR = [AVAR tau(avar_min)];
    
    x_hat = irreg_moving_avg_estimator(t, timestamps(recent), y(recent), tau, x_hat);
    [~, mse_min] = min(abs(x_hat - signal_ref(k).*ones(1,numel(tau))));
    MSE = [MSE tau(mse_min)];
    
    opt_est = [opt_est x_hat(mse_min)];
    t_est = [t_est t];
    avar_est = [avar_est x_hat(avar_min)];
    
    if mod(k,1) == 0
    subplot(4,1,1)    
        scatter(timestamps(recent),y(recent),'.')
            xlabel('Time [s]')
            ylabel('Measurement $\theta$')
            xlim([0,1])
            ylim([y_min,y_max])
            hold on
        scatter(timestamps,y,'.','MarkerEdgeAlpha',.05,'MarkerEdgeColor','k')
            xlabel('Time [s]')
            ylabel('Measurement $\theta$')
            xlim([0,1])
            ylim([y_min,y_max])
            xline(t,'Color','k','LineWidth',2);
            xline(t-tau(avar_min),'Color',c(1,:),'LineWidth',2);
            xline(t-tau(mse_min),'Color',c(2,:),'LineWidth',2);
            grid on
            hold off
            box on
    
    subplot(4,1,2)
        plot(tau,avar,'LineWidth',1.5)
%         hold on
%         plot(tau,avar)
        xlabel('Window length $\tau [s]$')
        ylabel('AVAR $\sigma^2_A$')
        xlim([0,0.5])
        ylim([0,200])
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        grid on
        box on
        
       
        subplot(4,1,3)
        t_prime = linspace(horizon,t,numel(AVAR));
        plot(t_prime,AVAR,'Color',c(1,:),'LineWidth',1.5)
        hold on
        plot(t_prime,movmean(MSE,[1 0]),'Color',c(2,:),'LineWidth',1.5)
        grid on
        xlim([0,1])
        ylim([-0.005,0.05])
        xlabel('Time [s]')
        ylabel('Window length $\tau [s]$')
        legend('$\tau_{avar}$','$\tau_{opt}$')
%         set(gca,'yscale','log')
        
        hold off
        pause(0.001)
        
        subplot(4,1,4)
            hold on
            plot(t_est,avar_est,'Color',c(1,:),'LineWidth',1.5);
            plot(t_est,opt_est,'Color',c(2,:)','LineWidth',1.5);
%             plot(t_est, signal_ref(1:ll),'Color','k','LineWidth',1.5);
            
            xlabel('Time [s]')
            ylabel('Estimated Signal')
            legend('avar est.','opt est.')
            xlim([0,1])
            hold off
            grid on
         
            k
    end
end



