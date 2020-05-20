function [avar] = AVAR(x,tau)
n = numel(x);
avar = [];
counter = 0;
for i=1:numel(tau)
    T = tau(i);
    Th = floor(T/2);
    avar_sum = 0;
    
    for m= 1+Th:n-3*Th
        avar_sum = avar_sum + (mean(x(m-Th:m+Th))-mean(x(m+Th:m+3*Th)))^2;
        counter = counter +1;
    end
    avar = [avar 0.5*avar_sum./counter];
end



