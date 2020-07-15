function y = OU(n)

t_start = 0;          %simulation start time
t_end = 1;          %simuation end time
dt = 1/n;            %time step
tau = 1/n * 2;            %relaxation time
c = 1;                %diffusion constant
x0 = 0;               %initial value for stochastic variable x
mu = 0;               %mean of stochatic process x
y0 = 0;               %initial value for integral x 
start_dist = -2.0;    %start of OU pdf 
end_dist = 2.0;       %end of OU pdf
%time
T = t_start:dt:t_end;
%compute x and y
i = 1;
x(1) = x0; 
y(1) = y0;
for t=t_start+dt:dt:t_end
   i = i + 1; 
   r1 = randn;
   r2 = randn;
   x(i) = x(i-1)*exp(-dt/tau) + sqrt((c*tau*0.5)*(1-(exp(-dt/tau))^2))*r1;
   y(i) = y(i-1) + x(i-1)*tau*(1-exp(-dt/tau))+sqrt((c*tau^3*(dt/tau-2*(1-exp(-dt/tau))+0.5*(1-exp(-2*dt/tau))))-((0.5*c*tau^2)*(1-exp(-dt/tau))^2)^2/((c*tau/2)*(1-exp(-2*dt/tau))))*r2+((0.5*c*tau^2)*(1-exp(-dt/tau))^2)/(sqrt((c*tau/2)*(1-(exp(-dt/tau))^2)))*r1;
end
%pdf for OU process
k = 0; j = start_dist:dt:end_dist;
for l=start_dist:dt:end_dist
    k = k + 1;
    p(k) = sqrt((1/tau)/(pi*c))*exp(-(1/tau)*(l-mu)^2/(c)); 
end
% plots
% figure;
% subplot(3,1,1)
% plot(T,x,'k-')
% xlabel('time')
% ylabel('x')
% subplot(3,1,2)
% plot(T,y,'k-')
% xlabel('time')
% ylabel('y')
% subplot(3,1,3)
% hold on
% histnorm(x,60)
% plot(j,p,'r-')
% xlabel('x')
% ylabel('probability')
% hold off