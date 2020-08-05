close all
clc
clear all


set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',16);

c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];

load('jerath_signal.mat');
n = 10000;
H = 100;
x_r = simout(1:n);

t = linspace(0,2*pi,n);
x_r = 1*sin(5*t)+2.0000*t+0.5000*sin(13*t);

for i=1:n
    if t(i) > 1 && t(i) < 2 || t(i) > 4 && t(i) < 5 
        x_r(i) = x_r(i) + 1*sin(50*t(i)) + 6;
    end
end

% x = zeros(n,1);

M = 20;

gamma = 1.2;

% for i= 1:M-1
%     tau = [floor(tau(1)/gamma) tau];
% end

c =1;



m = 50;
for i= 1:M-1
    if floor(m(1)/gamma) < 1
        break
    else
    m = [floor(m(1)/gamma) m];
    end
end

m = (1:50);
min_hist = [];
for k= 1:n
    
    sigma = 0.8;
    if c == 0
        x = [x normrnd(x_r(k),sigma)];
    else
        c = 0;
        x = normrnd(x_r(k),sigma);
    end
    

    
    avar = AVAR(x(max(1,end-H):end),m);
    
    
    
    [~,idx_avar] = min(avar);

    min_hist = [min_hist idx_avar];
    
%     subplot(2,1,1)
%     loglog(m,avar,'LineWidth',1.5)
% %     xlim([1 100])
%     xlabel('Window length [\# of samples]')
%     ylabel('Allan Variance')
% %     ylim([1e-11,1e-8]);
%     xlim([1 m(end)])
%     
%     subplot(2,1,2)
%     plot(x_r(1:k),'Color',c1,'LineWidth',2)
%     hold on
%     scatter((1:k),x,'MarkerEdgeColor',c2,'LineWidth',1.5)
%     legend('actual','measurement')
% %     ylim([-0.0003,0.0008])
%     xlabel('Time $[k]$')
%     ylabel('Parameter $\theta$')
%     title(['Noise scale: $\sigma=' num2str(sigma) '$'])
%     grid on
%     pause(0.05)
    
end


figure('Name','Signal + noise','Position',[20 50 800 500])
    scatter((1:k),x,'.','MarkerEdgeColor',c2,'LineWidth',1.5)
    hold on
    plot(x_r(1:k),'Color',c1,'LineWidth',2)
    legend('actual','measurement')
    xlabel('Time step $[k]$')
    ylabel('Parameter $\theta$')
    title(['Noise scale: $\sigma=' num2str(sigma) '$'])
    grid on

figure('Name','Window length vs. time','Position',[20 600 800 300])
    % plot(min_hist)
    % plot_with_bounds((1:n),min_hist,50,[0.9290, 0.6940, 0.1250]);
    plot(round(movmean(min_hist,100)),'LineWidth',2)
    xlabel('Time step [$k$]')
    ylabel('Window length $m$')
    grid on
