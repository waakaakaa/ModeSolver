function [ F ] = calBeta5LayerSymmetric( beta )
global lambda0;

global n1;
global n2;
global n3;

global a;
global d;

global m;
k0 = 2*pi/lambda0;

gama1 = sqrt(k0^2*n1^2 - beta.^2);
gama2 = sqrt(abs(beta.^2 - k0^2*n2^2));
gama3 = sqrt(beta.^2 - k0^2*n3^2);

T2 = gama2/gama1;
T3 = gama3/gama1;

delta2 = (T3 + T2*tanh(gama2*d)) / (T3*tanh(gama2*d) + T2) .* (beta>k0*n2)...
    +(T3 - T2*tan(gama2*d)) / (T3*tan(gama2*d) + T2) .* (beta<=k0*n2);
F = 2*gama1*a - m*pi - 2*atan(T2*delta2);
end