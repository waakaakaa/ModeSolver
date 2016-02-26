clc
clear

global lambda0;
global n0;
global refractiveIndex;
global nLast;
global thickness;
global mm;
global nn;
global neff_ridge;
global neff_slab;
global ridgeWidth


lambda0 = 1.55e-6;
n0 = 3.1562;
nLast = 1;
mm = 0;
nn = 0;
ridgeWidth = 3e-6;

figure(66);

%% ridge
refractiveIndex = [3.24896 3.24276 3.2703 ...
    3.4156 ...
    3.2703 3.24453 3.16853 3.45303  3.16724];
thickness = [0.01 0.06 0.06  ...
    0.0875 ...
    0.06 0.06 0.06 0.01  1.5]*1e-6;

% options = optimoptions('fsolve','Display','iter');
options = optimset('Display','none');
[neff_ridge,fval1] = fsolve(@calNeffQW_Ex_Step1, 3.2, options);

x = linspace(-2e-6, 4e-6, 100000);
E = plotExMultiLayer(neff_ridge ,x);
N = plotNMultiLayer(x);

subplot(1,2,1);hold on;box on;
plotyy(x*1e6,E,x*1e6,N);
title('ridge Ex and N');
hold off;

%% slab
refractiveIndex = [3.24896 3.24276 3.2703 ...
    3.4156 ...
    3.2703 3.24453 3.16853 3.45303];
thickness = [0.01 0.06 0.06  ...
    0.0875 ...
    0.06 0.06 0.06 0.01]*1e-6;

% options = optimoptions('fsolve','Display','iter');
[neff_slab,fval2] = fsolve(@calNeffQW_Ex_Step1, 3.16, options);

x = linspace(-2e-6, 4e-6, 100000);
E = plotExMultiLayer(neff_slab ,x);
N = plotNMultiLayer(x);

subplot(1,2,2);hold on;box on;
plotyy(x*1e6,E,x*1e6,N);
title('slab Ex and N');
hold off;

%% EIM
% options = optimoptions('fsolve','Display','iter');
[neff,fval3] = fsolve(@calNeffQW_Ex_Step2, 3.17, options);