function x_hat = reg_moving_avg_estimator(y, m)
n = numel(y);
x_hat = zeros(1,n);
for i= 1:n
    if i <= m
        x_hat(i) = mean(y(1:i));
    else
        x_hat(i) = mean(y(i-m+1:i));
    end
end
end
