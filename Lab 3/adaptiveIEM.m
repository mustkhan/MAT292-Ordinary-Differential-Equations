function [y,t] = adaptiveIEM(f, t0,tN,y0,h)

i=2;
h_init = h;

t_cur = t0;
t(1) = t0;
y(1) = y0;

tol = 1e-8;

while t_cur < tN
    f_full_step = 0.5*(f(t_cur,y(i-1)) + f(t_cur + h, y(i-1) + h*f(t_cur,y(i-1))));
    f_half_step = 0.5*(f(t_cur,y(i-1)) + f(t_cur + h/2, y(i-1) + (h/2)*f(t_cur,y(i-1))));
    
    Y = y(i-1) + f_full_step * h;
    Z_half = y(i-1) + f_half_step * h/2;
    
    f_half_step_2 = 0.5*(f(t_cur + h/2,Z_half) + f(t_cur + h, Z_half + (h/2)*f(t_cur + h/2,Z_half)));
    
    Z = Z_half + f_half_step_2 * h/2;
    
    D = Z-Y; % error estimate
    
    if abs(D) < tol
        y(i) = Z - D;   % local error O(h^3)
        t_cur = t_cur + h;
        t(i) = t_cur;
        i = i+1;
    end
    h = 0.9*h*min(max(tol/abs(D),0.3),2); 
end