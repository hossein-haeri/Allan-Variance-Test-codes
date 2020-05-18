function [avar] = AVAR(x,tau)
n = numel(x);
avar = [];
for i=1:numel(tau)
    T = tau(i);
    avar_sum = 0;
    for m= 1+T:n-2*T
        avar_sum = avar_sum + 0.5*(mean(x(m:m+2*T))-mean(x(m-T:m+T)))^2;
    end
    avar = [avar avar_sum/(n-2*T)];
end

