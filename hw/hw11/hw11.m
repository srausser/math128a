%% HW 11
% Sam Rausser
% 4/14/15
% 23485911

%% 2c)

fun = @(y, t) (y^2 + y)/t;
[w, t] = AdamsAdaptPC(fun, [1 3], -2, 10^-4, [.01 .4]);

%% 2a)
clear;
%y1 = y(1) - y(2) + 2;
%y2 = -y(1) + y(2) +4*t;

u1 = @(t) -0.5*exp(2*t) + t.^2 + 2*t - 0.5;
u2 = @(t) 0.5*exp(2*t) + t.^2 - 0.5;
exact = [u1(0:0.1:1); u2(0:0.1:1)]
[w, t] = rk4_systems(0, 1, 10, [-1; 0])
diff = abs(w - exact)

%% 4a)
clear;
%y''-3y'+2y = 6e^(-t)
% x1 = x(2);
% x2 = 6*exp(-t) + 3*x(2) - 2*x(1);
r = @(t) 2*exp(2*t) - exp(t) + exp(-t);
exact = r(0:0.1:1)
[w, t] = rk4_systems2(0, 1, 10, [2; 2])
diff = abs(w(1) - exact)

%% 6a)
clear;
u1 = @(t) -0.5*exp(2*t) + t.^2 + 2*t - 0.5;
u2 = @(t) 0.5*exp(2*t) + t.^2 - 0.5;
exact = [u1(0:0.1:1); u2(0:0.1:1)]
[w, t] = AdamsAdaptPC_systems(0, 1, 10, [-1; 0])
diff = abs(w - exact)