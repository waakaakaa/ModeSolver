function [ F ] = calNeffRidgeWG_Ey_Step3( neff )
global lambda0;
global N1;
global N2;
global a;
global m;
k0 = 2*pi/lambda0;
F = k0*sqrt(N1^2 - neff^2)*a...
    - m*pi...
    - 2*atan( sqrt(neff^2 - N2^2) / sqrt(N1^2 - neff^2) );

end