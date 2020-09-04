function x_hat = reg_moving_avg_estimator(t, timestamps , y, tau)
n = numel(t);
x_hat = zeros(1,n);
for i= 1:n
   window = y(timestamps<=t(i) & timestamps>t(i)-tau);
   x_hat(i)
    
end

end

