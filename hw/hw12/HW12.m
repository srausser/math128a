%% HW 12
% Sam Rausser
% 23485911

%% 2b)
clear t y;
a=0; b=1; N=10; y0=0; h=(b-a)/N;
t(1)=a;
y(1)=exp(1);
exact = @(t) exp(-10*t + 1) + t;
f = @(t, y) -10*y + 10*t + 1;
for n=1:N 
t(n+1)=t(n)+h;
y(n+1)=y(n)+h*f(t(n),y(n)); 
end
t = t';
approx = y(:); true = exact(a:h:b)'; diff = abs(y(:) - true);
table(t, approx, true, diff)

%% 4b)
clear t y;
f = @(t, y) -10*y + 10*t + 1;
[w, t_i] = RK4(f, [0 1], exp(1), 10);
diff = abs(w - true);
table(t_i, w, true, diff)

%% 6b)
clear t y w;
f = @(t, y) -10*y + 10*t + 1;
[w, t_i] = Adams4PC(f, [0 1], exp(1), 10);
diff = abs(w - true);
table(t_i, w, true, diff)

%% 8b)
clear t y w;
df = @(t, y) -10;
[w, t_i, flg] = TrapNewton(f, df, [0 1], exp(1), 10, 100, 10^-5);
diff = abs(w - true);
table(t_i, w, true, diff)
