function [ F ] = calNeffChannelWG_Ey_Step1( neff )
global lambda0;
global n1;
global n2;
global n4;
global b;
global n;
k0 = 2*pi/lambda0;
F = k0*sqrt(n1^2 - neff^2)*b...
    - n*pi...
    - atan( (n1^2*sqrt(neff^2 - n2^2)) / (n2^2*sqrt(n1^2 - neff^2)) )...
    - atan( (n1^2*sqrt(neff^2 - n4^2)) / (n4^2*sqrt(n1^2 - neff^2)) );

end