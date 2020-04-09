clc
close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

window_length_list = [0.01 0.1 1.0 10];
mycolor =  lines(7);
colors_p = parula(numel(window_length_list)+1);
colors_p = colors_p(1:end-1,:);

dt = 0.01;             % [sec]

for w= 1:numel(window_length_list)
x = [.5 0]';           % Initial state
ref = [0 0]';           % Reference state
b = 0.2;



window_length = window_length_list(w);
X = [];
u_hist = [];
win = b;
u = 0;
count_cnt = 0;
count_meas = 0;
num_window_samples = floor(window_length/dt);


for t= 0:dt:10
%% update the states
   b = 8.0 + 5.5*sin(2*t);
%    b = 0.5;
   x_dot = dyn_model(x,u,b);
   x = x + x_dot.*dt;

%    ref = [0.3*cos(t) 0]';
   if t > 2, ref = [0.5 0]'; end
   if t > pi, ref = [-0.3*sin(t) 0]'; end
   if t > 3*pi, ref = [0.0 0]'; end
   
%% make the noisy measurements

%    b_meas = b + 10.9*(rand()-0.5);
    b_meas = b + normrnd(0, 5);
  

%% take moving average of the state as the estimated state
   if size(win,2) < num_window_samples
       win = [win b_meas];
   else
       win = circshift(win,-1);
       win(:,end) = b_meas;
   end
   b_est = mean(win,2);
   
   X = [X [x; b; b_meas; b_est; ref]];
%% linearize the system with the estimated parameter at the current operating point
   [A, B] = make_linear_at(x, u, b_est);
   
%% obtain the input signal from controller
   u = controller(A, B, x, ref);



end

%% visualize the results
figure(1)
subplot(2,2,w)
hold on
t = dt*(1:1:size(X,2));
% p1 = plot(t, X(1,:),'DisplayName','$\theta$','LineWidth',1.5);
% p2 = plot(t, X(2,:),'DisplayName','$\dot{\theta}$','LineWidth',1.5);
p4 = plot(t, X(4,:),'DisplayName','$b_{meas}$','LineWidth',1.2, 'Color', mycolor(3,:),'LineStyle','-');
% p4 = scatter(t, X(4,:), 'Marker','.', 'DisplayName','$b_{meas}$', 'MarkerFaceColor', mycolor(4,:));
p5 = plot(t, X(5,:),'DisplayName','$b_{est}$','LineWidth',1.2, 'Color', mycolor(1,:));
p3 = plot(t, X(3,:),'DisplayName','$b_{true}$','LineWidth',1.2, 'Color', mycolor(2,:), 'LineStyle','--');
% p1.Color(4) = 0.5;
% p2.Color(4) = 0.5;
p3.Color(4) = 0.9;
p4.Color(4) = 0.7;
p5.Color(4) = 0.9;
hl = legend('show');
set(hl, 'Interpreter','latex','FontSize',12)
xlabel('$t [s]$','Interpreter','latex','FontSize',12)
ylabel('friction coefficient $b$','Interpreter','latex','FontSize',12)
title(['Window length = ' num2str(window_length) 's' ' (' num2str(num_window_samples) ' records' ')'],'Interpreter','latex','FontSize',14)
grid on




figure(2)
subplot(1,2,1)
hold on
t = dt*(1:1:size(X,2));
p1 = plot(t, X(1,:),'LineWidth',1.5,'DisplayName','$\theta$', 'Color',colors_p(w,:));
if w == numel(window_length_list)
pt1 = plot(t, X(6,:),'LineWidth',1.5,'DisplayName','$ref$', 'Color',mycolor(2,:), 'LineStyle','--');
end
p1.Color(4) = 0.8;


subplot(1,2,2)
hold on
t = dt*(1:1:size(X,2));
p2 = plot(t, X(2,:),'LineWidth',1.5,'DisplayName','$\dot{\theta}$', 'Color',colors_p(w,:));
if w == numel(window_length_list)
pt2 = plot(t, X(7,:),'LineWidth',1.5,'DisplayName','$ref$', 'Color',mycolor(2,:), 'LineStyle','--');
end
p2.Color(4) = 0.8;






end



figure(2)
subplot(1,2,1)
xlabel('$t$ [s]','Interpreter','latex','FontSize',14)
ylabel('$\theta$ [rad]','Interpreter','latex','FontSize',14)
grid on

subplot(1,2,2)
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

function u = controller(A, B, x, ref)
U = 200; % controller limit
    p = [-10 -8];
    K = place(A, B, p);
    u = -K*(x-ref);
    u = min(max(-U,u),U);
end



function x_dot = dyn_model(x, u, b)
g = 9.81;
m = 1;
l = 1;
I = m*l^2;
x_dot = [0 0]';
x_dot(1) = x(2);
x_dot(2) = (u+m*g*sin(x(1)))/I  - b*x(2);
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


