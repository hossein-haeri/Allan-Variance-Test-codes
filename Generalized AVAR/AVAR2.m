function [avar] = AVAR2(t,x,tau)
n = numel(x);
avar = [];
counter = 0;
t_min = min(t);
t_max = max(t);

t_range = t_max-t_min;

for i=1:numel(tau)
    T = tau(i);
    
    avar_sum = 0;
    
    while j < n
        t(j)
    end
    avar = [avar 0.5*avar_sum./counter];
end



