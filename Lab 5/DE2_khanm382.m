function y = DE2_khanm382(p,q,g,t0,tN,y0,y1,h)
t = t0:h:tN;
N = length(t);
y = zeros(1, N); %Solutions
y(1) = y0; %Initial condition
y(2) = y0 + y1*h; %Euler step
for i = 3:N-1
    y(i) = (h*h)*(g(t(i-1))-q(t(i-1))*y(i-1)-p(t(i-1))*((y(i-1)-y(i-2))/h)) + 2*y(i-1) - y(i-2);
end
