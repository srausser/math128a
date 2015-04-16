function f = vdpode_test(t,y)
f = [ y(2); sin(t) - y(1) - (y(1)^2 - 1)*y(2) ];
end