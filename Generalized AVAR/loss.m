
function e = loss(actual,estimate,N)
    e = zeros(1,N);
    for n= 1:N
%         e(n) = immse(actual,estimate(n,:));
        e(n) = sum(abs(actual-estimate(n,:)))/numel(actual);
    end

end



