%      n4
% n3   n1   n5
%      n2

% y
% |
% |
% |
% -------x

%      n4
% 
%      n1  (height = b)
% 
%      n2

%     (width = a)
% 
% n3       N1       n5

clc
clear

global lambda0
global n1
global n2
global n3
global n4
global n5
global m
global n
global a
global b

n1 = 1.4949;
n2 = 1.465;
n3 = 1.465;
n4 = 1.465;
n5 = 1.465;
a = 6e-6;
b = 6e-6;
lambda0 = 1.55e-6;
m = 2;
n = 1;

k0 = 2*pi/lambda0;
% options = optimoptions('fsolve','Display','iter');
options = optimset('Display','none');

[kx,fval1] = fsolve(@calKx_Ey, 10e5, options);
[ky,fval2] = fsolve(@calKy_Ey, 10e5, options);

kz = sqrt(k0^2*n1^2 - kx^2 - ky^2);

neff = kz/k0

x = linspace(-10e-6, 10e-6, 100);
y = linspace(-10e-6, 10e-6, 100);
[ Ey_x, Ey_y ] = plotEy(x, y, kx,ky);
figure(111);hold on;box on;
plot(x, Ey_x,'r-')
plot(y, Ey_y,'b*')
xlim([-10e-6 10e-6])
hold off;

Ey = Ey_x'*Ey_y;
figure(222);hold on;box on;
pcolor(Ey);shading interp;colorbar;
hold off;