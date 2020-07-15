close all
clc
clear all

figure('Position',[20 20 800 800])
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',16);

c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];

load('jerath_signal.mat');
n = 1000;
x_r = simout(1:n);

t = linspace(0,2*pi,1000);
x_r = 0.0002*sin(3*t)+0.0000*t+0.0000*sin(7*t);
% x = zeros(n,1);

M = 20;

gamma = 1.2;

% for i= 1:M-1
%     tau = [floor(tau(1)/gamma) tau];
% end

c =1;

for k= 1:n
    
    sigma = 0.0002*sin(k/300);
    if c == 0
        x = [x normrnd(x_r(k),sigma)];
    else
        c = 0;
        x = normrnd(x_r(k),sigma);
    end
    
    tau = floor(numel(x)/1.5);
    for i= 1:M-1
        if floor(tau(1)/gamma) < 1
            break
        else
        tau = [floor(tau(1)/gamma) tau];
        end
    end
    
    avar = AVAR(x(max(1,end-floor(n/10)):end),tau);
    
    subplot(2,1,1)
    loglog(tau,avar,'LineWidth',1.5)
%     xlim([1 100])
    xlabel('Window length [\# of samples]')
    ylabel('Allan Variance')
    ylim([1e-11,1e-8]);
%     xlim([1 13])
    
    subplot(2,1,2)
    plot(x_r(1:k),'Color',c1,'LineWidth',2)
    hold on
    scatter((1:k),x,'MarkerEdgeColor',c2,'LineWidth',1.5)
    legend('actual','measurement')
    ylim([-0.0003,0.0008])
    xlabel('Time $[k]$')
    ylabel('Parameter $\theta$')
    title(['Noise scale: $\sigma=' num2str(sigma) '$'])
    pause(0.05)
    
end




