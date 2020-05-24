close all
clear all
clc
set(groot,'defaulttextinterpreter','tex');
set(groot, 'defaultAxesFontSize',16);
N = 1;
color = lines(N+1);
num_samples = 100;
t = linspace(1,num_samples,num_samples);
% t = sort(rand(1,num_samples));


q = floor(numel(t)/2);
y = 0.05*t.^1.8;

data = [];
b = [20 10];
for n= 1:N
%    data = [data; y + b(n)*(rand(1,numel(t))-0.5)];
    data = [data; y + b(n).*(normrnd(0,1,[1 numel(t)]))];
end


hold on

c = lines(4);
% for i= 1:num_samples
%     xline(t(i),'-','Color',[.8 .8 .8],'LineWidth',1);
% %     alpha(0.4)
% end
% 


% for i= 1:num_samples
%     xline(t(i),'-','Color',[.8 .8 .8],'LineWidth',1);
% %     alpha(0.4)
% end
% 

for i= 1:num_samples/10
    k = i*10;
%     xline(t(k),'-','Color',[.8 .8 .8],'LineWidth',1);
    
%     scatter(t(k)-5,mean(data(k-9:k)),50,'O','MarkerFaceColor',c(2,:),'MarkerEdgeColor',c(2,:));
end



% scatter(t,data,10,'filled','MarkerFaceAlpha',0.8)

plot(t,data,'LineWidth',2)


set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

xlabel('Time')
ylabel('x')
