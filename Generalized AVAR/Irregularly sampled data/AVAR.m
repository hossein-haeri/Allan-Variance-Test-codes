function [avar] = AVAR(x,tau)
n = numel(x);
avar = [];
for i=1:numel(tau)
    T = tau(i);
    avar_sum = 0;
    c = 0;
    for k= 1:n-2*T
        avar_sum = avar_sum + 0.5*(mean(x(k:k+T))-mean(x(k+T:k+2*T)))^2;
        c = c + 1;
    end
    avar = [avar avar_sum./c];
end

