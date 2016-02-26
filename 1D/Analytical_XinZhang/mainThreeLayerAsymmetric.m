clc;
clear;

global lambda0;
lambda0 = 0.6328e-6;
k0 = 2*pi/lambda0;

global n1;
n1 = 1.568;

global n2;
n2 = 1.538;

global n3;
n3 = 1;

global b;
b = 4e-6;

global m;

count = 1;
neff = zeros(1,4);
for m = [0 1 2 3]
    %     options = optimoptions('fsolve','Display','none');
    options = optimset('Display','none');
    [x,fval] = fsolve(@calBeta3LayerAsymmetric, k0*(n1-0.01),options);

    beta0 = real(x(1));
    gama1 = sqrt(k0^2*n1^2 - beta0^2);
    gama2 = sqrt(beta0^2 - k0^2*n2^2);
    gama3 = sqrt(beta0^2 - k0^2*n3^2);
    neff(1, count) = beta0/k0;
    
    x = linspace(-b,2*b,100);
    
    Ey = exp(gama2*x)                                                  .* (x < 0)+...
        (cos(gama1*x) + (gama2/gama1)*sin(gama1*x))                    .* (x >=0 & x <= b)+...
        (cos(gama1*b) + (gama2/gama1)*sin(gama1*b))*exp(-gama3*(x-b))  .* (x > b);
    
    subplot(2,2,count);
    plot(x*1e6,Ey);
    count = count + 1;
end
figure(2);plot(neff, 'ro--');