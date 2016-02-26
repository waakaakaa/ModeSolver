clc
clear

global lambda0;
global n0;
global refractiveIndex;
global nLast;
global thickness;
global m;

%% 3 layers
% lambda0 = 0.6328e-6;
% n0 = 1.538;
% refractiveIndex = [1.598];
% nLast = 1.538;
% thickness = [4]*1e-6;
% m = 0;
%
% options = optimoptions('fsolve','Display','iter');
% [neff,fval] = fsolve(@calNeffMultiLayer, max(refractiveIndex)-0.05, options);
%
% % neff0 = 1.5665
%
% x = linspace(-thickness, 2*thickness, 1000);
% E = plotExMultiLayer( neff ,x);
% plot(x*1e6,E,'r.-')

%% 4 layers
% lambda0 = 0.98e-6;
% n0 = 3.5;
% refractiveIndex = [3.55 3.45];
% nLast = 1;
% thickness = [1 2]*1e-6;
% m = 0;
%
% options = optimoptions('fsolve','Display','iter');
% [neff,fval] = fsolve(@calNeffMultiLayer, max(refractiveIndex)-0.025, options);
%
% x = linspace(-2e-6, 4e-6, 100);
% E = plotExMultiLayer( neff ,x);
% plot(x*1e6,E,'r.-')

%% 5 layers
% lambda0 = 1.55e-6;
% n0 = 3.5214;
% refractiveIndex = [3.5941 3.9262 3.5941];
% nLast = 3.5214;
% thickness = [0.5 0.5 0.5]*1e-6;
% m = 0;
%
% options = optimoptions('fsolve','Display','iter');
% [neff,fval] = fsolve(@calNeffMultiLayer, max(refractiveIndex)-0.2, options);
%
% x = linspace(-2e-6, 4e-6, 1000);
% E = plotExMultiLayer( neff ,x);
% plot(x*1e6,E,'r.-')

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
[neff,fval] = fsolve(@calNeffMultiLayer, 3.2, options);

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
[neff,fval] = fsolve(@calNeffMultiLayer, 3.2, options);

E = plotExMultiLayer( neff ,x);
plot(x*1e6,E,'b-')
hold off;