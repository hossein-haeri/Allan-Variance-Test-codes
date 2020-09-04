
function x_hat = irreg_moving_avg_estimator(t, timestamps, y, tau_list, x_hat_prev)

x_hat = zeros(1,numel(tau_list));
% x_hat_prev = 0;
for i= 1:numel(tau_list)
    tau = tau_list(i);
    window = y(timestamps<=t & timestamps>t-tau);
    if numel(window) == 0
        x_hat(i) = x_hat_prev(i);
    else
        x_hat(i) = mean(window);
    end
end
end

