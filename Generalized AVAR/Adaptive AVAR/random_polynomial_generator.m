function [y, x_r] = random_polynomial_generator(t_r, t, order)


c = 10*(rand(1,order+1)-0.5);

y = c(end);
x_r = c(end);
for r= 1:order
    y = y + c(r).*t.^r;
    x_r = x_r + c(r).*t_r.^r;
end

% y = y(1:end-1);
end