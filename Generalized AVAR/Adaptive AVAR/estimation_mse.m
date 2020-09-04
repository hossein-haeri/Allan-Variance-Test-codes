function mse = estimation_mse(t_r, x_r, t_hat, x_hat)

n = numel(x_r);
x = zeros(1,n);
j = 1;
% resample estimated signal up to the original signal
for i= 1:n
   if t_r(i) >= t_hat(j) && t_r(i) < t_hat(j+1)
      x(i) = x_hat(j);
   elseif t_r(i) < t_hat(j)
      x(i) = initial_guess;
   else 
      x(i) = x_hat(j+1);
      j = j + 1;
   end
end

mse = mean((x_r - x).^2);

end

