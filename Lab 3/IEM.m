function y = IEM(g,t0,tN,y0,h) %this is my IEM solver
t = t0:h:tN;
N = length(t);
y = zeros(1, N);
y(1) = y0;
for i = 1:N-1
    k1 = g(t(i), y(i));
    k2 = g(t(i)+h, y(i)+(h*k1));
    y(i+1) = y(i) + h*((k1+k2)/2);
end
    
    
    