function [a,b] = lin_reg(x,y)

%  [a,b] = lin_reg(x,y)
%  y = a + bx
x = x(:);
y = y(:);
Sx = sum(x);
Sy = sum(y);
Sxx = sum(x.*x);
Sxy = sum(x.*y);
N = length(x);
Mx = Sx/N;
My = Sy/N;

b = (Sx*My - Sxy)/(Sx*Mx - Sxx);
a = My - b*Mx;
