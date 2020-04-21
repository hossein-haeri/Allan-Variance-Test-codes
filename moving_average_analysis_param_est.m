clc
close all
clear all

%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
window_length_list = [0.01 0.05 0.1 0.2 1.0 5];
dt = 0.01;             % [sec]
x0 = [.2 0]';          % Initial state



set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','tex');
set(groot, 'defaultLegendInterpreter','latex');
mycolor =  lines(7);
colors_p = parula(numel(window_length_list)+1);
colors_p = colors_p(1:end-1,:);



for w= 1:numel(window_length_list)
x = x0;           % set x to the initial state


window_length = window_length_list(w);
X = [];
u_hist = [];
win = [];
E = [];
horizon = 50;
u = 0;
count_cnt = 0;
count_meas = 0;
num_window_samples = floor(window_length/dt);
E_norm = [];

for t= 0:dt:50
%% update the parameter value and states
   theta = 2.0 + 1*sin(2*t);
%    b = 2.0;
   x_dot = dyn_model(x,u,theta);
   x = x + x_dot.*dt;
   
   
%% build a reference trajectory
    if t >= 0, ref = [0 0]'; end
%     if t > 5,  ref = [0.2 0]'; end
%     if t > 2*pi, ref = [-0.1*sin(2*t)+0.1 0]'; end
  
%% make a noisy measurements
%    theta_meas = theta + 10.0*(rand()-0.5);
    theta_meas = theta + normrnd(0, 1.5);
  

%% take moving average of the state as an estimate of the state
   if size(win,2) < num_window_samples
       win = [win theta_meas];
   else
       win = circshift(win,-1);
       win(:,end) = theta_meas;
   end
   theta_est = mean(win,2);
   
   
%% linearize the system with the estimated parameter at the current operating point
   [A_est, B_est] = make_linear_at(x, u, theta_est);
   [A, B] = make_linear_at(x, u, theta);
   
%% obtain the input signal from controller
   [u, K] = controller(A_est, B_est, x, ref);

    A_cl = A - B*K;
    A_cl_est = A_est - B_est*K;

    lam = eig(A_cl);
    lam_est = eig(A_cl_est);
    e = norm(lam - [-1; 3])^2;
    
   if numel(E) < horizon
       E = [E e];
   else
       E = circshift(E,-1);
       E(:,end) = e;
   end
   eigen_variance = sum(E)/horizon^2;
    
    
%% record the data into X    
    X = [X [x; theta; theta_meas; theta_est; ref; u; eigen_variance]];
end

%% visualize the results
figure(1)
subplot(2,3,w)
hold on
t = dt*(1:1:size(X,2));
p1 = plot(t, X(4,:),'DisplayName','$b_{meas}$','LineWidth',1.2, 'Color', mycolor(3,:),'LineStyle','-');
p2 = plot(t, X(5,:),'DisplayName','$b_{est}$','LineWidth',1.2, 'Color', mycolor(1,:));
p3 = plot(t, X(3,:),'DisplayName','$b_{true}$','LineWidth',1.2, 'Color', mycolor(2,:), 'LineStyle','--');
p1.Color(4) = 0.7;
p2.Color(4) = 0.9;
p3.Color(4) = 0.9;
hl = legend('show');
set(hl, 'Interpreter','latex','FontSize',12)
xlabel('$t [s]$','Interpreter','latex','FontSize',12)
ylabel('friction coefficient $b$','Interpreter','latex','FontSize',12)
title(['Window length = ' num2str(window_length) 's' ' (' num2str(num_window_samples) ' records' ')'],'Interpreter','latex','FontSize',12)
grid on




figure(2)
subplot(2,1,1)
hold on
t = dt*(1:1:size(X,2));
p1 = plot(t, X(1,:),'LineWidth',1.0,'DisplayName','$\theta$', 'Color',colors_p(w,:));
if w == numel(window_length_list)
pt1 = plot(t, X(6,:),'LineWidth',1.0,'DisplayName','$ref$', 'Color',mycolor(2,:), 'LineStyle','--');
end
p1.Color(4) = 0.8;
title('Pendulum Position','Interpreter','latex','FontSize',12)


subplot(2,1,2)
hold on
t = dt*(1:1:size(X,2));
p2 = plot(t, X(2,:),'LineWidth',1.0,'DisplayName','$\dot{\theta}$', 'Color',colors_p(w,:));
if w == numel(window_length_list)
pt2 = plot(t, X(7,:),'LineWidth',1.0,'DisplayName','$ref$', 'Color',mycolor(2,:), 'LineStyle','--');
end
p2.Color(4) = 0.8;
title('Pendulum Velocity','Interpreter','latex','FontSize',12)



figure(3)
subplot(2,1,1)
hold on
e_norm = sqrt((X(1,:)-X(6,:)).^2+(X(2,:)-X(7,:)).^2);
E_norm = cumsum(e_norm);
pe = plot(t, E_norm,'LineWidth',1.0,'DisplayName','$\theta$', 'Color',colors_p(w,:));
pe.Color(4) = 0.8;

figure(3)
subplot(2,1,2)
hold on
pu = plot(t, X(8,:),'LineWidth',1.0,'DisplayName','$\theta$', 'Color',colors_p(w,:));
pu.Color(4) = 0.8;



figure(4)
pv = plot_with_bounds(t,X(9,:),30, colors_p(w,:));
pv.Color(4) = 0.8;
hold on

end



figure(2)

    subplot(2,1,1)
        xlabel('$t$ [s]','Interpreter','latex','FontSize',14)
        ylabel('$\theta$ [rad]','Interpreter','latex','FontSize',14)
        grid on
        colormap(colors_p);
        Tk = linspace(1/(2*numel(window_length_list)),1 - 1/(2*numel(window_length_list)),numel(window_length_list));
        cb = colorbar('Ticks',  Tk, 'TickLabels',string(window_length_list));
        cb.Label.String = 'window length [s]';
        cb.FontSize = 12;
        cb.Label.Interpreter = 'latex';
        cb.TickLabelInterpreter = 'latex';

    subplot(2,1,2)
        xlabel('$t$ [s]','Interpreter','latex','FontSize',14)
        ylabel('$\dot{\theta}$ [rad/s]','Interpreter','latex','FontSize',14)
        grid on
        colormap(colors_p);
        Tk = linspace(1/(2*numel(window_length_list)),1 - 1/(2*numel(window_length_list)),numel(window_length_list));
        cb = colorbar('Ticks',  Tk, 'TickLabels',string(window_length_list));
        cb.Label.String = 'window length [s]';
        cb.FontSize = 12;
        cb.Label.Interpreter = 'latex';
        cb.TickLabelInterpreter = 'latex';


figure(3)
    subplot(2,1,1)
        xlabel('$t$ [s]','Interpreter','latex','FontSize',14)
        ylabel('cumulative $|e|$','Interpreter','latex','FontSize',14)
        colormap(colors_p);
        Tk = linspace(1/(2*numel(window_length_list)),1 - 1/(2*numel(window_length_list)),numel(window_length_list));
        cb = colorbar('Ticks',  Tk, 'TickLabels',string(window_length_list));
        cb.Label.String = 'window length [s]';
        cb.FontSize = 12;
        cb.Label.Interpreter = 'latex';
        cb.TickLabelInterpreter = 'latex';
        grid on


figure(3)
    subplot(2,1,2)
        xlabel('$t$ [s]','Interpreter','latex','FontSize',14)
        ylabel('input signa $u$ [Nm] ','Interpreter','latex','FontSize',14)
        colormap(colors_p);
        Tk = linspace(1/(2*numel(window_length_list)),1 - 1/(2*numel(window_length_list)),numel(window_length_list));
        cb = colorbar('Ticks', Tk, 'TickLabels',string(window_length_list));
        cb.Label.String = 'window length [s]';
        cb.FontSize = 12;
        cb.Label.Interpreter = 'latex';
        cb.TickLabelInterpreter = 'latex';
        grid on


figure(4)
    xlabel('$t$ [s]','Interpreter','latex')
    ylabel('Eigenvalues variance','Interpreter','latex')
    colormap(colors_p);
    Tk = linspace(1/(2*numel(window_length_list)),1 - 1/(2*numel(window_length_list)),numel(window_length_list));
    cb = colorbar('Ticks', Tk, 'TickLabels',string(window_length_list));
    cb.Label.String = 'window length [s]';
    cb.FontSize = 12;
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';
    grid on

function [u,K] = controller(A, B, x, ref)
    U = 50; % controller limit
    p = [-1 -3];
    K = place(A, B, p);
%     K = [19 10];
    u = -K*(x-ref);
    u = min(max(-U,u),U);
end



function x_dot = dyn_model(x, u, m)
g = 9.81;
b = 0.2;
% m = 1;
l = 1;
I = m*l^2;
x_dot = [0 0]';
x_dot(1) = x(2);
x_dot(2) = (u+m*g*sin(x(1)))/I  - (b*x(2)).^3;
end




function [A, B] = make_linear_at(x0, u0, b)
n = numel(x0);
m = numel(u0);
A = zeros(n);
B = zeros(n,m);
ep = 0.01;
for i= 1:n
    ep_vector = zeros(n,1);
    ep_vector(i) = ep;
    A(:,i) = (dyn_model(x0+ep_vector, u0, b) - dyn_model(x0, u0, b))./ep;
end
for i= 1:m
    ep_vector = zeros(m,1);
    ep_vector(i) = ep;
    B(:,i) = (dyn_model(x0, u0+ep_vector, b) - dyn_model(x0, u0, b))./ep;
end
end


