function [rts,info] = cubic23485911(C)
%project 1 
%   Want to find roots of cubic polynomial
%   https://math.berkeley.edu/~mgu/MA128A2015S/Prog1.pdf
% func = poly2sym(C);
% disp(func);

a = C(1);
b = C(2);
c = C(3);
d = C(4);
TOL = 10^(-6);
N0 = 10^6;
p0 = -2;
p1 = 2;
p2 = 0;
    
if (a == 0)
    if (b == 0)
        if (c == 0)
            info = sprintf('no roots as polynomial entered is constant on R');
        end
        rts(1) = -d/c;
        info = sprintf('one root located at %d', rts(1));
    end
    a = C(2);
    b = C(3);
    c = C(4);
    rts(1) = (-b + sqrt(b^2 - 4*a*c))/(2*a);
    rts(2) = (-b - sqrt(b^2 - 4*a*c))/(2*a);
    info = sprintf('two roots located at %d and %d', rts(1), rts(2));
    return;
end
    
% step1:    
h1 = p1 - p0; 
h2 = p2 - p1;
delta1 = (polyval(C, p1) - polyval(C, p0))/h1; 
delta2 = (polyval(C, p2) - polyval(C, p1))/h2; 
d = (delta2 - delta1)/(h2 + h1);
i = 3;

% step 2 
while (i <= N0) % do Steps 3-7
	% step 3
    b = delta2 + h2 * d;
    D = (b^2 -4*polyval(C, p2)*d)^(1/2);

	% step 4	
    if (abs(b - D) < abs(b + D))
        E = b + D;
    else
        E = b - D;
    end

	% step 5
    h = -2*polyval(C, p2)/E;
	p = p2 + h;
    
	% step 6
    if (abs(h) < TOL)
		rts(1) = p;
        break;
    end

	% step 7 (prepare for next iteration.)
    p0 = p1; 	
	p1 = p2; 
    p2 = p;
	h1 = p1 - p0;
	h2 = p2 - p1;
	delta1 = (polyval(C, p1) - polyval(C, p0))/h1;
	delta2 = (polyval(C, p2) - polyval(C, p1))/h2;
	d = (delta1 - delta2)/(h2 + h1);
	i = i + 1;
end

b = [1 rts(1)];
[q, r] = deconv(C,b);
disp(r);
a = q(1);
b = q(2);
c = q(3);
rts(2) = (-b + sqrt(b^2 - 4*a*c))/(2*a);
rts(3) = (-b - sqrt(b^2 - 4*a*c))/(2*a);
info = sprintf('three roots located at %d, %d, and %d', rts(1), rts(2), rts(3));
return;
end

