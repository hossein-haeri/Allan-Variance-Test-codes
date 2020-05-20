close all
clear all
clc

N = 4;
color = lines(N+1);
num_samples = 2000;
t = linspace(0,10,num_samples);
q = floor(numel(t)/2);
y = 80*sin(t).^1;

data = [];
sigma = [50 100 120 30];
for n= 1:N
%    data = [data; y + b(n)*(rand(1,numel(t))-0.5)];
    data = [data; y + (normrnd(0,sigma(n),[1 numel(t)]))];
end

W = 1./sigma.^2;
data = [data; W*data./sum(W)];


N = N+1;


y = y(floor(numel(t)/2):end);



M = 20;
gamma = 1.2;
tau = floor(numel(t)/gamma^5);

while tau(1) > 5
    tau = [floor(tau(1)/gamma) tau];
end

avar = zeros(N,numel(tau));
for n= 1:N
avar(n,:) = AVAR(data(n,:),tau);
end


%% Simulation
L = [];
I = findmin(avar,N);

m_opt = zeros(n,1);
for n= 1:N
    m_opt(n) = tau(I(n));
end
W = fliplr(m_opt')/norm(m_opt,1);

optimal_y_hat = zeros(N,q+1);

for j= 1:numel(tau)
    y_hat_hist = [];
    m = tau(j);

    for k= q:numel(t)
        y_hat = zeros(N,1);
        for n= 1:N
            y_hat(n) = mean(data(n,max(1,k-m):k));

        end
            y_hat_hist = [y_hat_hist y_hat];
            
    end

    L = [L; loss(y,y_hat_hist,N)];
    
    for n= 1:N
        if j == I(n)
            optimal_y_hat(n,:) = y_hat_hist(n,:);
        end
    end
end



% agg_y_hat = mean(W*optimal_y_hat,1);



t = t(q:end);    
data = data(:,q:end);

figure(1)
    hold on
    plot(t,y,'LineWidth',2,'Color','k','LineStyle','--','DisplayName','actual')
    alpha(0.2)
    for n= 1:N
        if n == N
            plot(t,optimal_y_hat(n,:),'LineWidth',1.0,'Color',color(n,:),'LineStyle','-','DisplayName','Fusion Estimation')
        else
            plot(t,optimal_y_hat(n,:),'LineWidth',1.0,'Color',color(n,:),'LineStyle','-','DisplayName',['Estimation ' int2str(n)])
        end
        scatter(t,data(n,:),'.','MarkerEdgeColor',color(n,:),'MarkerEdgeAlpha',0.2,'DisplayName','measurement','HandleVisibility','off')
    end
    xlabel('time [s]')
    ylabel('\theta')
    legend
    
%     plot(t,agg_y_hat,'LineWidth',1,'Color',color(N+1,:),'LineStyle','-','DisplayName','aggrigated estimation')
    xlim([min(t) max(t)])

figure(2)
subplot(2,1,1)
    hold on
    for n= 1:N
        if n==N
            plot(tau,avar(n,:),'Color',color(n,:),'LineWidth',1.5,'DisplayName','Fusion')
        else
            plot(tau,avar(n,:),'Color',color(n,:),'LineWidth',1.5,'DisplayName',['Sensor ' int2str(n)])
        end
        scatter(tau(I(n)),avar(n,I(n)),100,'*','MarkerEdgeColor',color(n,:),'HandleVisibility','off')
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Window length m [samples]')
    ylabel('Allan variance  \sigma_\theta^2')
    legend
    grid on
    
subplot(2,1,2)
hold on
%     plot(tau,repmat(loss(y,agg_y_hat,1),[1 numel(tau)]),'LineWidth',2,'Color',color(N+1,:),'LineStyle','--')
    for n= 1:N
        if n == N
        plot(tau,L(:,n),'LineWidth',1.5,'Color',color(n,:),'DisplayName','Fusion')
        else
        plot(tau,L(:,n),'LineWidth',1.5,'Color',color(n,:),'DisplayName',['Sensor ' int2str(n)])
        end
        [v,idx]=min(L(:,n));
        scatter(tau(idx),v,100,'*','MarkerEdgeColor',color(n,:),'HandleVisibility','off')
    end
    
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Window length m [samples]')
    ylabel('MSE')
    grid on
    
function I = findmin(avar,N)
    I = zeros(1,N);
    d = diff(avar,1,2);
    for n= 1:N
        I(n) = find(diff(d(n,:)>=0),1)+1;
    end
end


