function [w,t] = RK4(FunFcnIn, Intv, alpha, N)
% On input: 
%   FunFcnIn is the name of function to be integrated
%   interv is the interval to be integrated over
% 
%   The problem: y' = f(t,y), y(a) = alpha, a<= t <= b
%   where Intv = [a b], and a call to FunFcnIn with 
%   argument (t, y) returns f(t,y).
%
%   RK4 uses the explicit 4th order Runge-Kutta method 
%   to approximate y at N+1 equi-spaced points on [a, b].
%
%
% On output
%   w contains the approximations at the mesh points.
%
% Written by Ming Gu for Math 128A, Fall 2008
% 
[FunFcn,msg] = fcnchk(FunFcnIn,0);
if ~isempty(msg)
    error('InvalidFUN',msg);
end
a    = Intv(1);
b    = Intv(2);
h    = (b-a)/N;
w    = zeros(N+1,1);
w(1) = alpha;
t    = a + (0:N)'*h;
h2   = h/2;
for i = 1:N
    k1 = h* FunFcn(t(i),w(i));
    k2 = h* FunFcn(t(i)+h2,w(i)+k1/2);
    k3 = h* FunFcn(t(i)+h2,w(i)+k2/2);
    k4 = h* FunFcn(t(i)+h,w(i)+k3);
    w(i+1) = w(i) + (k1+2*k2+2*k3+k4)/6;
end