function [w,t,flg] = RKFv_test(Intv, alpha, tol, stepsize)
% On input: 
%   FunFcnIn is the name of function to be integrated
%   interv is the interval to be integrated over
% 
%   The problem: y' = f(t,y), y(a) = alpha, a<= t <= b
%   where Intv = [a b], and a call to FunFcnIn with 
%   argument (t, y) returns f(t,y).
%
%   RKF uses the Runge-Kutta-Fehlberg method to solve
%   the above problem, to a given tolerance. 
%   hmin and hmax are the minimum and maximum step sizes.
% 
% this is the vector version of RKF.
%
% On output
%   t contains the (unequi-spaced) mesh points, w the
%   function values at these points. 
%   flg is success flag. If flg == 0 method succeeded,
%   and if flg ~=0 method failed.
%
% Written by Ming Gu for Math 128A, Fall 2008
% 
% [FunFcn,msg] = fcnchk(FunFcnIn,0);
% if ~isempty(msg)
%     error('InvalidFUN',msg);
% end
flg  = 0;
a    = Intv(1);
b    = Intv(2);
hmin = stepsize(1);
hmax = stepsize(2);
h    = min(hmax,b-a);
if (h<=0)
    msg = ['Illegal hmax or interval'];
    error('InvalidFUN',msg);
end
w    = alpha;
t    = a;
c    = [0;      1/4; 3/8;       12/13;      1;  1/2];
d    = [25/216  0    1408/2565  2197/4104 -1/5  0   ];
r    = [1/360   0   -128/4275  -2197/75240 1/50 2/55];
KC   = [zeros(1,5); 
        1/4,       zeros(1,4); 
        3/32       9/32       zeros(1,3);
        1932/2197 -7200/2197, 7296/2197, 0,         0;
        439/216   -8          3680/513  -845/4104,  0;
          -8/27    2         -3544/2565  1859/4104 -11/40];
K    = zeros(length(alpha),6);
%
% main loop
%
while (flg == 0)
    ti   = t(end);
    wi   = w(:,end);
    for j = 1:6
        K(:,j) = h*vdpode_test(ti+c(j)*h,wi+K(:,1:5)*KC(j,:)');
    end
%
% accept approximation
%
    R = max(eps,norm(K*r')/h);
    if (R <= tol)
        t = [t; ti+h];
        w = [w, wi+K*d'];
    end
%
% reset stepsize.
%
    delta = 0.84*(tol/R)^(1/4);
    delta = min(4,max(0.1,delta));
    h     = min(hmax, delta * h);
    h     = min(h, b-t(end));
    if (abs(h)<eps)
        return;
    end
    if (h < hmin)
        flg = 1;
        return;
    end
end
