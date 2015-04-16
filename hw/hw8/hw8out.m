%% HW #8
%  Sam Rausser
%  23485911
%  3/17/15

%% 4.6
%% 4d)
fun = @(x) exp(x).*(6*cos(4*x)+4*sin(6*x))
[Int,flg,fcnt,level] = AdaptSimpson(fun, [0, pi/2], 0.00001, 100)
clear

%% 5b)
fun = @(x) x.*sin(x.^2);
q = integral(fun, 0, pi);
n = 5;
while(abs(q - CompSimpson(fun, n, [0 , pi])) > 0.000001)
    n = n + 2;
end
n
[Int,flg,fcnt,level] = AdaptSimpson(fun, [0 , pi], 0.000001, 10)
clear

%% 6b)

x = 0.1:pi/100:2;
y = cos(1./x);
% figure
% plot(x, y)
y = @(x) cos(1./x);
[Int,flg,fcnt,level] = AdaptSimpson(y, [0.1, 2], .001, 10)
clear

%% 7)
m = 1;
k = 9;
F_0 = 1;
w = 2;
w_0 = sqrt(k/m);
% t = 0:pi/20:2*pi
% y = F_0/(m*(w_0^2-w^2)).*(cos(w*t)-cos(w_0*t))
% plot(t, y)
u = @(t) F_0/(m*(w_0^2-w^2)).*(cos(w*t)-cos(w_0*t));
[Int,flg,fcnt,level] = AdaptSimpson(u, [0, 2*pi], .0001, 10)
clear


%% 8a)
m = 1;
k = 9;
F_0 = 1;
c = 10;
w = 2;
w_0 = sqrt(k/m);
r_1 = (-c+sqrt(c^2-4*w_0^2*m^2))/2*m;
r_2 = (-c-sqrt(c^2-4*w_0^2*m^2))/2*m;
syms c_1 c_2 t
u =  c_1*exp(r_1*t) + c_2*exp(r_2*t) + (F_0/(c^2*w^2+m^2*(w_0^2-w^2)^2)).*(c*w*sin(w*t)+m*(w_0^2-w^2)*cos(w*t))
diff_u = diff(u)
u_0 = vpa(subs(u, t, 0))
diff_u_0 = vpa(subs(diff_u, t, 0))
% c_1 = solve(u_0 == 0) 
c_1 = -1/40;
c_2 = 9/680;
u = @(t) c_1*exp(r_1*t) + c_2*exp(r_2*t) + (F_0/(c^2*w^2+m^2*(w_0^2-w^2)^2)).*(c*w*sin(w*t)+m*(w_0^2-w^2)*cos(w*t))
[Int,flg,fcnt,level] = AdaptSimpson(u, [0, 2*pi], 0.0001, 10)
clear

%% 4.7
%% 1b)
a = 0;
b = 1;
syms t x
f = @(x) (x.^2).*exp(-1*x);
actual = int(f, x, a, b)
x = (t+1)/2;
% dx = 1/2 * dt;
[c,x_val] = Legendre(2)
ans = vpa((1/2) * (c(1) * subs(f(x), t, x_val(1)) + c(2) * subs(f(x), t, x_val(2))))
diff = vpa(abs(actual - ans))
clear

%% 2b)
a = 0;
b = 1;
syms t x
f = @(x) (x.^2).*exp(-1*x);
actual = int(f, x, a, b)
x = (t+1)/2;
% dx = 1/2 * dt;
[c,x_val] = Legendre(3)
ans = vpa((1/2) * (c(1) * subs(f(x), t, x_val(1)) + c(2) * subs(f(x), t, x_val(2)) + c(3) * subs(f(x), t, x_val(3))))
diff = vpa(abs(actual - ans))
clear





