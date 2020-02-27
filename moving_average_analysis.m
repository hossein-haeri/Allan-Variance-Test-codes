clc
close all



x = [1.5 0]';           % Initial state
ref = [0 0]';           % Reference state

controller_freq = 10;   % [Hz]
sensor_freq = 500;      % [Hz]
dt = 0.001;             % [sec]
window_length = 0.5;    % [sec]


X = [];
u_hist = [];
win = x;
u = 0;
x_meas = x;
count_cnt = 0;
count_meas = 0;

for t= 0:dt:10
   %% update the states
   x_dot = dyn_model(x,u);
   x = x + x_dot.*dt;
   X = [X x];
   
   %% make the noisy measurements
   count_meas = count_meas + 1;
   if count_meas >= floor(1/(sensor_freq*dt))
        count_meas = 0;
       x_meas(1) = x(1) + 0.5*(rand()-0.5);
       x_meas(2) = x(2) + 0.1*(rand()-0.5);
   end
   %% take moving average of the state as the estimated state
   if size(win,2) < floor(window_length/dt)
       win = [win x_meas];
   else
       win = circshift(win,-1,2);
       win(:,end) = x_meas;
   end
   x_est = mean(win,2);
   %% linearize the system and get the control signal u
   [A, B] = make_linear_at(x_est, u);
   count_cnt = count_cnt + 1;
   if count_cnt >= floor(1/(controller_freq*dt))
    count_cnt = 0;
       u = controller(A, B, x_est, ref);
   end
%     u_hist = [u_hist u];

end

%% visualize the results
hold on
t = dt*(1:1:size(X,2));
% plot(t, u_hist,'LineWidth',0.1);
p1 = plot(t, X(1,:),'DisplayName','$\theta$','LineWidth',1.5);
p2 = plot(t, X(2,:),'DisplayName','$\dot{\theta}$','LineWidth',1.5);
p2.Color(4) = 0.5;
p1.Color(4) = 0.8;
hl = legend('show');
set(hl, 'Interpreter','latex','FontSize',12)
xlabel('t [s]')
ylabel('states')
title(['Window length = ' num2str(window_length) ' s'],'FontSize',10)
grid on




function u = controller(A, B, x, ref)
U = 10;
%     p = [-15-6i -15+6i];
    p = [-10 -5];
    K = place(A, B, p);
    u = -K*(x-ref);
    u = min(max(-U,u),U);
end

function x_dot = dyn_model(x, u)
g = 9.81;
m = 1;
b = 0.2;
l = 1;
I = m*l^2;
x_dot = [0 0]';
x_dot(1) = x(2);
x_dot(2) = (u+m*g*sin(x(1)))/I  - b*x(2);
end


function [A, B] = make_linear_at(x0, u0)
n = numel(x0);
m = numel(u0);
A = zeros(n);
B = zeros(n,m);
ep = 0.01;
for i= 1:n
    ep_vector = zeros(n,1);
    ep_vector(i) = ep;
    A(:,i) = (dyn_model(x0+ep_vector, u0) - dyn_model(x0, u0))./ep;
end
for i= 1:m
    ep_vector = zeros(m,1);
    ep_vector(i) = ep;
    B(:,i) = (dyn_model(x0, u0+ep_vector) - dyn_model(x0, u0))./ep;
end
end



