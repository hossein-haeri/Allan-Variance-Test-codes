function y = random_walk(size)
y = zeros(size,1);
for i= 2:size
    a = rand;
%     if a >= 0.5
%     y(i) = y(i-1) + 0.5;
%     else
%     y(i) = y(i-1) - 0.5;
%     end
    y(i) = y(i-1) + a-0.5;
end

end