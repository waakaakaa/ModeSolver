%      n3
%      n1   
%      n2

% y
% |
% |
% |
% -------x

%     (width = a)
% 
% height = 
% b-h      b        b-h
%
% N2       N1       N2

clc
clear

global lambda0
global n1
global n2
global n3
global m
global n
global a
global b
global h
global N1
global N2

n1 = 1.4949;
n2 = 1.465;
n3 = 1.465;
lambda0 = 1.55e-6;
m = 0;
n = 0;

bm = (0.1:0.1:6)*1e-6;
neff = zeros(1,length(bm));
count = 1;
figure(111);hold on;box on;
for b = bm
    a = 1.5*b;
    h = 0.25*b;
    
    %     options = optimoptions('fsolve','Display','iter');
    options = optimset('Display','none');
    [N1,fval1] = fsolve(@calNeffRidgeWG_Ey_Step1, 1.466, options);
    [N2,fval2] = fsolve(@calNeffRidgeWG_Ey_Step2, 1.466, options);
    [N, fval3] = fsolve(@calNeffRidgeWG_Ey_Step3, 1.466, options);
    neff(1,count) = N;
    count = count + 1;
end
plot(bm*1e6, neff,'r.-')
xlim([0 6])
ylim([1.465 1.491])