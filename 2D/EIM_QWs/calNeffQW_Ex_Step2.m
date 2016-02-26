function [ F ] = calNeffQW_Ex_Step2( neff )
global lambda0;
global neff_ridge;
global neff_slab;
global ridgeWidth;
global mm;
k0 = 2*pi/lambda0;
F = k0*sqrt(neff_ridge^2 - neff^2)*ridgeWidth...
    - mm*pi...
    - 2*atan( (neff_ridge^2*sqrt(neff^2 - neff_slab^2)) ...
    / (neff_slab^2*sqrt(neff_ridge^2 - neff^2)) );

end