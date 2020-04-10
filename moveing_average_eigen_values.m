clear all, clc

clc
close all
fig1 = figure(1);
a1 = axes;
fig2 = figure(2);
a2 = axes;
% window_length = 0.1;    % [sec]

window_length_list = [0.01 0.05 0.1 0.2 1.0 10];
% window_length_list = [0.005];

[A, B] = make_linear_at([0 0]', 0);
for w= 1:numel(window_length_list)

window_length = window_length_list(w);
x = [0.2 0]';           % Initial state
ref = [0 0]';           % Reference state

controller_freq = 10;   % [Hz]
sensor_freq = 200;      % [Hz]
dt = 0.01;             % [sec]



                         

X = [];
u_hist = [];
win = x;
u = 0;
K = [0 0];
x_meas = x;
count_cnt = 0;
count_meas = 0;
num_window_samples = floor(window_length/dt);
A_bar = zeros(2*num_window_samples);
B_bar = zeros(2*num_window_samples,1);
L_bar = zeros(1,2*num_window_samples);
K_bar = zeros(1,2*num_window_samples);
EV = [];
% K_bar = zeros(


    for t= 0:dt:5
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
       if size(win,2) < num_window_samples
           win = [win x_meas];
       else
           win = circshift(win,-1,2);
           win(:,end) = x_meas;
       end
       x_est = mean(win,2);




       %% linearize the system and get the control signal u
       [A, B] = make_linear_at(x_est, u);

       % get sensor dynamics
    %    c = moving_avg_tf(dt,num_window_samples);
    %    H = tf2ss(c{1},c{2});

       count_cnt = count_cnt + 1;
       if count_cnt >= floor(1/(controller_freq*dt))
        count_cnt = 0;
           [u, K] = controller(A, B, x_est, ref);
       end
    %     u_hist = [u_hist u];   


       x_bar = reshape(win,numel(win),1);
       A_bar(1:2,1:2) = A;
       for i= 1:num_window_samples-1
            A_bar(2*i+1:2*i+2,2*i-1:2*i) = eye(2);
       end
       B_bar(1:2,1) = B;
       L_bar = eye(2*num_window_samples)/num_window_samples;
       K_bar(1:2) = K;
       A_cl = A_bar - B_bar*K_bar*L_bar;

       ev = eig(A_cl);
       EV = [EV ev(end-1:end)];

    end
label = char(['w=' num2str(window_length) 's' ' (' num2str(num_window_samples) ' record(s)' ')']);
t = dt*(1:1:size(X,2));
figure(1)
plot(a1,t, real(EV(2,:)),'DisplayName',label, 'LineWidth', 1)
figure(2)
plot(a2,t, real(EV(1,:)),'DisplayName',label, 'LineWidth', 1)

hold(a1,'on')
hold(a2,'on')

end
                      


figure(1)
grid on
a1 = legend('show');
set(a1, 'Interpreter','latex','FontSize',12)
a1 = xlabel('t [s]');
title(['Window length = ' num2str(window_length) ' s'],'FontSize',10)

figure(2)
grid on
a2 = legend('show');
set(a2, 'Interpreter','latex','FontSize',12)
a2 = xlabel('t [s]');
title(['Window length = ' num2str(window_length) ' s'],'FontSize',10)


% % visualize the results
% hold on
% t = dt*(1:1:size(X,2));
% plot(t, u_hist,'LineWidth',0.1);
% p1 = plot(t, X(1,:),'DisplayName','$\theta$','LineWidth',1.5);
% p2 = plot(t, X(2,:),'DisplayName','$\dot{\theta}$','LineWidth',1.5);
% p2.Color(4) = 0.5;
% p1.Color(4) = 0.8;
% hl = legend('show');
% set(hl, 'Interpreter','latex','FontSize',12)
% xlabel('t [s]')
% ylabel('states')
% title(['Window length = ' num2str(window_length) ' s'],'FontSize',10)
% grid on
% 

%%

function [u, K] = controller(A, B, x, ref)
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

function coeff = moving_avg_tf(T, N)
num_terms = 5;
Num = [];
Den = [];

    for i= 0:num_terms
        Num = [Num -T*N.^i/factorial(i)];
        Den = [Den -T.^i/factorial(i)];
    end
Num(1) = 1 - Num(1);
Den(1) = (1 - Den(1)).*N;
% mv_tf = tf(Num, Den);
coeff = {Num, Den};
end

