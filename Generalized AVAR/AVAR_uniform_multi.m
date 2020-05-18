close all
clear all
clc

N = 2;
color = lines(N+1);

t = linspace(0,10,1000);
q = floor(numel(t)/2);
y = 2*t.^2+1;

data = [];
b = [60 20];
for n= 1:N
%    data = [data; y + b(n)*(rand(1,numel(t))-0.5)];
    data = [data; y + b(n).*(normrnd(0,1,[1 numel(t)]))];
end


y = y(floor(numel(t)/2):end);



M = 20;
gamma = 1.5;
tau = floor(numel(t)/gamma^3);

while tau(1) > 1
    tau = [floor(tau(1)/gamma) tau];
end

avar = zeros(N,numel(tau));
for n= 1:N
avar(n,:) = AVAR(data(n,:),tau);
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
        plot(t,optimal_y_hat(n,:),'LineWidth',1,'Color',color(n,:),'LineStyle','-','DisplayName','Avar estimation')
        scatter(t,data(n,:),'.','MarkerEdgeColor',color(n,:),'MarkerEdgeAlpha',0.4,'DisplayName','measurement')
    end
    
%     plot(t,agg_y_hat,'LineWidth',1,'Color',color(N+1,:),'LineStyle','-','DisplayName','aggrigated estimation')


figure(2)
subplot(2,1,1)
    hold on
    for n= 1:N
        plot(tau,avar(n,:),'Color',color(n,:),'LineWidth',2)
        scatter(tau(I(n)),avar(n,I(n)),100,'*','MarkerEdgeColor',color(n,:))
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('window length (# of samples)')
    ylabel('Allan variance')
    grid
    
subplot(2,1,2)
hold on
%     plot(tau,repmat(loss(y,agg_y_hat,1),[1 numel(tau)]),'LineWidth',2,'Color',color(N+1,:),'LineStyle','--')
    for n= 1:N
        plot(tau,L(:,n),'LineWidth',2,'Color',color(n,:))
        [v,idx]=min(L(:,n));
        scatter(tau(idx),v,100,'*','MarkerEdgeColor',color(n,:))
    end
    
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('window length (# of samples)')
    ylabel('MSE')
    
    grid
    
function I = findmin(avar,N)
    I = zeros(1,N);
    d = diff(avar,1,2);
    for n= 1:N
        I(n) = find(diff(d(n,:)>=0),1)+1;
    end
end

