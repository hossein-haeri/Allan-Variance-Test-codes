clear all
close all
clc

X = [0];
u = 0;
for i=1:10
    x = 2*X(end)+1 + 2*u;
    X = [X x];
    
end

figure('Position',[50 50 400 300])
plot(X,'LineWidth',1.5)
grid on
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$x_k$','Interpreter','latex','FontSize',16)
