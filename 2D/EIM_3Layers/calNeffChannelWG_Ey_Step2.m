function [ F ] = calNeffChannelWG_Ey_Step2( neff )
global lambda0;
global N1;
global n3;
global n5;
global a;
global m;
k0 = 2*pi/lambda0;
F = k0*sqrt(N1^2 - neff^2)*a...
    - m*pi...
    - atan( sqrt(neff^2 - n3^2) / sqrt(N1^2 - neff^2) )...
    - atan( sqrt(neff^2 - n5^2) / sqrt(N1^2 - neff^2) );

end