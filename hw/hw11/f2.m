function dy = f2(t,y)
dy = zeros(2,1); % a column vector
dy(1) = y(2);
dy(2) = 6*exp(-t) + 3*y(2) - 2*y(1);
end