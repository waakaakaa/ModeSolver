function [ F ] = calBeta3LayerSymmetric( beta )
global lambda0;
global n1;
global n2;
global a;
global m;
k0 = 2*pi/lambda0;
gama1 = sqrt(k0^2*n1^2 - beta.^2);
gama2 = sqrt(beta.^2 - k0^2*n2^2);
F = 2*gama1*a - m*pi - 2*atan(gama2/gama1);
end