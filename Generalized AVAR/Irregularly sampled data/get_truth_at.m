function y = get_truth_at(t,tout, simout)
    
    y = zeros(1,numel(t));
    for i=1:numel(t)
        [~,indx] = min(abs(tout-t(i)));
        y(i) = simout(indx);
        % y = 1.0*t + 1*sin(t/1).^1+2;

    end
end