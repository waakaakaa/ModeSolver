clc;
clear;

global lambda0;
lambda0 = 0.98e-6;
k0 = 2*pi/lambda0;

global n1;
n1 = 3.55;

global n2;
n2 = 3.5;

global n3;
n3 = 3.45;

global n4;
n4 = 1;

global b1;
b1m = (0.3: 0.1: 1)*1e-6;
% b1 = 1e-6;

global b3;
b3m = (0: 0.1: 2.0)*1e-6;

global m;
m = 0;

%%
% countCol = 1;
% countRow = 1;
% neff = zeros(length(b1m), length(b3m));
%
% figure(1);
% hold on;
% box on;
% for b1 = b1m
%     for b3 = b3m
%         options = optimoptions('fsolve','Display','none');
%         [x,fval] = fsolve(@calBeta4LayerAsymmetric, k0*(n1-0.01), options);
%
%         beta0 = real(x(1));
%         neff(countCol, countRow) = beta0/k0;
%
%         countRow = countRow + 1;
%     end
%     plot(b3m*1e6, neff(countCol,:), 'r.-');
%     countCol = countCol + 1;
%     countRow = 1;
% end
% xlim([0 2]);
% ylim([3.5 3.54]);
% hold off;

% %%
% b1m = (0.3: 0.1: 3)*1e-6;
% neff = zeros(1, length(b1m));
% count = 1;
% plotType = 'rgb';
% countType = 1;
% figure(2);
% hold on;
% box on;
% for m = [0 1 2 3]
%     for b3 = [0 0.2 1.0]*1e-6
%         for b1 = b1m
%             options = optimoptions('fsolve','Display','none');
%             [x,fval] = fsolve(@calBeta4LayerAsymmetric, k0*(n1-0.01), options);
%
%             beta0 = real(x(1));
%             neff(1, count) = beta0/k0;
%             count = count + 1;
%         end
%         plot(b1m*1e6, neff(1,:), plotType(countType));
%         count = 1;
%         countType = countType + 1;
%     end
%     countType = 1;
% end
% xlim([0 3]);
% ylim([3.5 3.55]);
% hold off;


%%
% b3 = 0.2e-6;
% figure(3);hold on;box on;
% m = 0;
% for b1 = (0.5:0.1:1.5)*1e-6
%     options = optimoptions('fsolve','Display','iter');
%     [beta,fval] = fsolve(@calBeta4LayerAsymmetric, k0*(n1-0.001), options);
%     gama1 = sqrt(k0^2*n1^2 - beta.^2);
%     gama2 = sqrt(beta.^2 - k0^2*n2^2);
%     gama3 = sqrt(beta.^2 - k0^2*n3^2);
%     gama4 = sqrt(beta.^2 - k0^2*n4^2);
%     T2 = gama2/gama1;
%     T3 = gama3/gama1;
%     T4 = gama4/gama1;
%     A1 = 1;
%     B1 = T2;
%     A3 = A1*cos(gama1*b1)+B1*sin(gama1*b1);
%     B3 = (1/T3)*(-A1*sin(gama1*b1)+B1*cos(gama1*b1));
%     A4 = A3*cosh(gama3*b3)+B3*sinh(gama3*b3);
%     
%     x = linspace(-2e-6, 3e-6, 100);
%     E = exp(gama2*x).*(x<0) + ...
%         (A1*cos(gama1*x) + B1*sin(gama1*x)).*(x>=0&x<b1)+...
%         (A3*cosh(gama3*(x - b1)) + B3*sinh(gama3*(x - b1))).*(x>=b1&x<b3+b1)+...
%         A4*exp(-gama4*(x - b1 - b3)).*(x>=b1+b3);
%     plot(x,E);
% end
% m = 1;
% for b1 = (1.1:0.1:1.5)*1e-6
%     options = optimoptions('fsolve','Display','iter');
%     [beta,fval] = fsolve(@calBeta4LayerAsymmetric, k0*(n1-0.001), options);
%     gama1 = sqrt(k0^2*n1^2 - beta.^2);
%     gama2 = sqrt(beta.^2 - k0^2*n2^2);
%     gama3 = sqrt(beta.^2 - k0^2*n3^2);
%     gama4 = sqrt(beta.^2 - k0^2*n4^2);
%     T2 = gama2/gama1;
%     T3 = gama3/gama1;
%     T4 = gama4/gama1;
%     A1 = 1;
%     B1 = T2;
%     A3 = A1*cos(gama1*b1)+B1*sin(gama1*b1);
%     B3 = (1/T3)*(-A1*sin(gama1*b1)+B1*cos(gama1*b1));
%     A4 = A3*cosh(gama3*b3)+B3*sinh(gama3*b3);
%     
%     x = linspace(-2e-6, 3e-6, 100);
%     E = exp(gama2*x).*(x<0) + ...
%         (A1*cos(gama1*x) + B1*sin(gama1*x)).*(x>=0&x<b1)+...
%         (A3*cosh(gama3*(x - b1)) + B3*sinh(gama3*(x - b1))).*(x>=b1&x<b3+b1)+...
%         A4*exp(-gama4*(x - b1 - b3)).*(x>=b1+b3);
%     plot(x,E);
% end
% hold off

%%
b1 = 1e-6;
b3 = 2e-6;
m = 0;
% options = optimoptions('fsolve','Display','iter');
options = optimset('Display','none');
[beta,fval] = fsolve(@calBeta4LayerAsymmetric, k0*(n1-0.001), options);
neff = beta/k0;