function [ Ey_x, Ey_y ] = plotEy(x, y, kx,ky)
global lambda0;
global n1;
global n2;
global n3;
global n4;
global n5;
global a;
global b;

k0 = 2*pi/lambda0;
k3x = sqrt(k0^2*(n1^2-n3^2) - kx^2);
k5x = sqrt(k0^2*(n1^2-n5^2) - kx^2);
k2y = sqrt(k0^2*(n1^2-n2^2) - ky^2);
k4y = sqrt(k0^2*(n1^2-n4^2) - ky^2);

kai = atan(k3x/kx) - kx*a/2;
Ey_x = cos(kx*a/2 + kai) .* exp(k3x*(x + a/2)) .* (x<-a/2)...
    + cos(kx*x - kai) .* (x<=a/2 & x>=-a/2)...
    + cos(kx*a/2 - kai) .* exp(-k5x*(x - a/2)) .* (x>a/2);

eta = ky*b/2 - atan((n1^2/n4^2)*(k4y/ky));
Ey_y = cos(ky*b/2 + eta) .* exp(k2y*(y + b/2)) .* (y<-b/2)...
    + cos(ky*x - eta) .* (y<=b/2 & y>=-b/2)...
    + cos(ky*b/2 - eta) .* exp(-k4y*(y - b/2)) .* (y>b/2);

end