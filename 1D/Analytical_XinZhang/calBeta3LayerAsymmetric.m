function [ F ] = calBeta3LayerAsymmetric( beta )
global lambda0;
global n1;
global n2;
global n3;
global b;
global m;
k0 = 2*pi/lambda0;
gama1 = sqrt(k0^2*n1^2 - beta.^2);
gama2 = sqrt(beta.^2 - k0^2*n2^2);
gama3 = sqrt(beta.^2 - k0^2*n3^2);
F = gama1*b - m*pi - atan(gama2./gama1) - atan(gama3./gama1);
end