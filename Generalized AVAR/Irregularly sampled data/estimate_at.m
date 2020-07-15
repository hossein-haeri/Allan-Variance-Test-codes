function [t_hat,y_hat] = estimate_at(data_t,data_y,tau,n)
T = data_t(end)-data_t(1);
dt = T/n;
y_hat = zeros(1,n);
t_hat = linspace(data_t(1),data_t(end),n);
t=0;
for i=1:n
    y_bar = mean(data_y(data_t<=t & data_t>t-tau));
    if ~isnan(y_bar)
        y_hat(i) = y_bar;
    elseif i>1
        y_hat(i) = y_hat(i-1);
    end
    t = t+dt;
end
end
