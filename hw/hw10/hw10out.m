%% HW #10
% Sam Rausser
% 23485911
% 4/7/15

%% 10a)
w = -1;
j = 0;
for t = 1.05:0.05:2
    w = w + 0.05*(1/(1 + 0.05*j)^2 - w./(1 + 0.05*j) - w^2 + 0.025*(-3/(1+0.05*j)^3 + (3*w^2)/(1+0.05*j) + 2*w^3));
    A(j+1,1) = j+1;
    A(j+1,2) = t;
    A(j+1,3) = w;
    A(j+1,4) = -1/t;
    A(j+1,5) = abs(A(j+1,3) - A(j+1, 4));
    j = j + 1;
end
T = array2table(A, 'VariableNames',{'i' 't_i' 'w_i', 'exact_val', 'error'})


