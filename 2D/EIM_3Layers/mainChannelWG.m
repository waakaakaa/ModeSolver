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
global N1

n1 = 1.5;
n2 = 1;
n3 = 1;
n4 = 1;
n5 = 1;
a = 2e-6;
b = 1e-6;
% lambda0 = 2e-6;

m = 1;
n = 0;

B = (0:0.1:4);
P = zeros(1,length(B));

count = 1;
figure(111);hold on;box on;
for m = [0 1 2]
    for lambda0 = 2*b*sqrt(n1^2 - n2^2)./B;
        %         options = optimoptions('fsolve','Display','iter');
        options = optimset('Display','none');
        [N1,fval1] = fsolve(@calNeffChannelWG_Ey_Step1, 1.01, options);
        [N,fval2] = fsolve(@calNeffChannelWG_Ey_Step2, 1.01, options);
        P(1,count) = (N^2 - n2^2)/(n1^2 - n2^2);
        count = count + 1;
    end
    plot(B,P,'r.-');
    count = 1;
end
xlim([0 4])
ylim([0 1])