function [y, x_r] = random_sine_generator(t_r, t, freq_list)


r = numel(freq_list);
c = 2*(rand(1,r)-0.5);

y = 0;
x_r = 0;
for i= 1:r
    y = y + c(r).*sin(freq_list(i).*t) + 2*t;
    x_r = x_r + c(r).*sin(freq_list(i).*t_r);
end

% plot(y);
end