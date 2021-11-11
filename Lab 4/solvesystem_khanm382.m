function [t, M] = solvesystem_khanm382(f,g,t0,tN,x0,h)
t = t0:h:tN;
N = length(t);
x1 = zeros(1, N);
x2 = zeros(1, N);
x1(1) = x0(1);
x2(1) = x0(2);
for i = 1:N-1
    k1 = f(t(i), x1(i), x2(i));
    u1 = g(t(i), x1(i), x2(i));
    
    k2 = f(t(i)+h, x1(i)+(h*k1), x2(i)+(h*u1));
    x1(i+1) = x1(i) + h*((k1+k2)/2);
    
    u2 = g(t(i)+h, x1(i)+(h*k1), x2(i)+(h*u1));
    x2(i+1) = x2(i) + h*((u1+u2)/2);
end
M = [x1; x2];
    