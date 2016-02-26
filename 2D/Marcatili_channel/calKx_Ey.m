function [ F ] = calKx_Ey( kx )
global lambda0;
global n1;
global n3;
global n5;
global a;
global m;
k0 = 2*pi/lambda0;
k3x = sqrt(k0^2*(n1^2-n3^2) - kx^2);
k5x = sqrt(k0^2*(n1^2-n5^2) - kx^2);
F = kx*a - m*pi - atan(k3x/kx) - atan(k5x/kx);
end