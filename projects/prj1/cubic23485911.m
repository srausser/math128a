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
TOL = 10^(-14);
N0 = 10^3;
% disp(C);
    
if (abs(a) == 0)
    if (abs(b) == 0)
        if (abs(c) == 0)
            info = sprintf('no roots as polynomial entered is constant on R');
            rts(1) = nan;
            return;
        end
        rts(1) = -d/c;
        info = sprintf('one root located at %d', rts(1));
        return;
    end
    a = C(2);
    b = C(3);
    c = C(4);
    rts(1) = (-b + sqrt(b^2 - 4*a*c))/(2*a);
    rts(2) = (-b - sqrt(b^2 - 4*a*c))/(2*a);
    info = sprintf('two roots located at %d and %d', rts(1), rts(2));
    return;
end

% der = polyder(C);
% ader = der(1);
% bder = der(2);
% cder = der(3);
% r1der = (-bder + sqrt(bder^2 - 4*ader*cder))/(2*ader);
% r2der = (-bder - sqrt(bder^2 - 4*ader*cder))/(2*ader);

% if (isreal(r1der) && isreal(r2der)) %if not true then only one real root
%     p0 = -1000;
%     p1 = 1000;
%     p2 = 0;
% elseif (a > 0)
%     p0 = -1000;
%     p1 = 100;
%     p2 = -450;
% else
%     p0 = -100;
%     p1 = 1000;
%     p2 = 450;
% end

randoms = randn(1,3);
% randoms = [-2, 0, 2];

while (polyval(C, randoms(2)) - polyval(C, randoms(1)) == 0 || polyval(C, randoms(3)) - polyval(C, randoms(2)) == 0 || polyval(C, randoms(3)) - polyval(C, randoms(1)) == 0)
    randoms = 2*randoms;
end

% while (polyval(C, randoms(3)) - polyval(C, randoms(2)) == 0)
%     randoms = 2*randoms;
% end
% 
% while (polyval(C, randoms(3)) - polyval(C, randoms(1)) == 0)
%     randoms = 2*randoms;
% end

% disp(randoms);

p0 = randoms(1);
p1 = randoms(2);
p2 = randoms(3);




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

disp(i);
if (not(exist('rts', 'var')))
%     disp('                          not exists');
    [rts,info] = cubic23485911(C);
    return;
end

b = [1 -rts(1)];
[q, r] = deconv(C,b);
% disp(r);
a = q(1);
b = q(2);
c = q(3);
rts(2) = (-b + sqrt(b^2 - 4*a*c))/(2*a);
rts(3) = (-b - sqrt(b^2 - 4*a*c))/(2*a);
info = sprintf('three roots located at %d, %d, and %d', rts(1), rts(2), rts(3));
return;
end

