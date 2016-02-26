function [ F ] = calKy_Ey( ky )
global lambda0;
global n1;
global n2;
global n4;
global b;
global n;
k0 = 2*pi/lambda0;
k2y = sqrt(k0^2*(n1^2-n2^2) - ky^2);
k4y = sqrt(k0^2*(n1^2-n4^2) - ky^2);
F = ky*b - n*pi - atan((n1^2/n2^2)*(k2y/ky)) - atan((n1^2/n4^2)*(k4y/ky));
end