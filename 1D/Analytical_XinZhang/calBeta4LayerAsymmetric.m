function [ F ] = calBeta4LayerAsymmetric( beta )
global lambda0;
global n1;
global n2;
global n3;
global n4;
global b1;
global b3;
global m;
k0 = 2*pi/lambda0;
gama1 = sqrt(k0^2*n1^2 - beta.^2);
gama2 = sqrt(beta.^2 - k0^2*n2^2);
gama3 = sqrt(beta.^2 - k0^2*n3^2);
gama4 = sqrt(beta.^2 - k0^2*n4^2);
T2 = gama2/gama1;
T3 = gama3/gama1;
T4 = gama4/gama1;
delta3 = ((T4+T3)+(T4-T3)*exp(-2*gama3*b3))/((T4+T3)-(T4-T3)*exp(-2*gama3*b3));
F = gama1*b1 - m*pi - atan(T2) - atan(T3*delta3);
end