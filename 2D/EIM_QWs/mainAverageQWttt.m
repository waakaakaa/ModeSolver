clc
clear

global lambda0;
global n0;
global refractiveIndex;
global nLast;
global thickness;
global m;

%% 5 QW approximation
lambda0 = 1.55e-6;
n0 = 3.1562;
refractiveIndex = [3.24896 3.24276 3.2703 ...
    3.4156 ...
    3.2703 3.24453 3.16853 3.45303  3.16724];
nLast = 1;
thickness = [0.01 0.06 0.06 ...
    0.0875 ...
    0.06 0.06 0.06 0.01  1.5]*1e-6;
m = 0;

% options = optimoptions('fsolve','Display','iter');
options = optimset('Display','none');
[neff,fval] = fsolve(@calNeffQW_Ex_Step1, 3.2, options);

figure(66);hold on;box on;
x = linspace(-2e-6, 4e-6, 100);
E = plotExMultiLayer( neff ,x);
plot(x*1e6,E,'r*')

refractiveIndex = [3.24896 3.24276 3.2703 ...
    3.3647 3.5241 3.3647 3.5241 3.3647 3.5241 3.3647 3.5241 3.3647 3.5241 3.3647 ...
    3.2703 3.24453 3.16853 3.45303  3.16724];
thickness = [0.01 0.06 0.06  ...
    0.01 0.0055 0.01 0.0055 0.01 0.0055 0.01 0.0055 0.01 0.0055 0.01   ...
    0.06 0.06 0.06 0.01  1.5]*1e-6;

% options = optimoptions('fsolve','Display','iter');
[neff2,fval2] = fsolve(@calNeffQW_Ex_Step2, 3.2, options);

E = plotExMultiLayer( neff ,x);
plot(x*1e6,E,'b-')
hold off;