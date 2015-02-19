function x0=bisection(func, a, b)
   
    fa=func(a);
    fb=func(b);
    if (fa*fb>0)
        error('Bisection requires a sign change on the initial interval.\n')
    end
    if (a>b)
        error('a should be less thanb.\n')
    end
    
    if (fa == 0)
        x0 = a;
        return;
    elseif (fb == 0)
        x0 = b;
        return;
    end
    
    % Assume fa and fb are both nonzero
    while (abs((a-b)/a) > eps)
        m = (a+b)/2;
        fm = func(m);
        if (fa*fm < 0)
            %that means [a, m] has a sign change. We should keep this
            %subinterval
            b = m;
            fb = fm;
        elseif (fb*fm < 0)
            %then [m, b] has a sign change
            a = m;
            fa = fm;
        else
            % We have the zero at m.
            x0 = m;
            return;
        end
    end
    x0 = (a+b)/2
end
