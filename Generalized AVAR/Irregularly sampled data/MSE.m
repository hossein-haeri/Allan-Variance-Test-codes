function L = MSE(t,y,tau_list,tout, simout)
L = [];
t_max = max(t);
t_min = min(t);
y_t = get_truth_at((t_min:0.01:t_max),tout, simout)';
initial_guess = 0;
% for each window length tau
for j= 1:numel(tau_list)
    y_hat_hist = [];
    tau = tau_list(j);
    % start the simulate
    for time=t_min:0.01:t_max
        % collect data points that fall into the window 
        valid_hist = y(t<=time & t>time-tau);
        % if window is not empty
        if ~isempty(valid_hist)
            % estimate the current value as the average over the window
            y_hat = mean(valid_hist);
        else
            % for empty window cases
            y_hat = initial_guess;
        end
        y_hat_hist = [y_hat_hist; y_hat];
    end
    % compute the MSE of the moving average estimate using window length tau
    l = sum((y_t-y_hat_hist).^2./numel(y_t));
    L = [L; l];
end
end

