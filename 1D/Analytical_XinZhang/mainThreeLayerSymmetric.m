clc;
clear;

global lambda0;
lambda0 = 0.6328e-6;
k0 = 2*pi/lambda0;

global n1;
n1 = 1.568;

global n2;
n2 = 1.538;

global a;
a = 2e-6;

global m;

count = 1;
neff = zeros(1,4);
for m = [0 1 2 3]
    phi = mod(m,2) * pi/2;
    
    beta = linspace(k0*n2,k0*n1,100);
    %     options = optimoptions('fsolve','Display','none');
    options = optimset('Display','none');
    [x,fval] = fsolve(@calBeta3LayerSymmetric,beta,options);
    
    beta0 = real(x(1));
    gama1 = sqrt(k0^2*n1^2 - beta0^2);
    gama2 = sqrt(beta0^2 - k0^2*n2^2);
    neff(1, count) = beta0/k0;
    
    x = linspace(-3*a,3*a,1000);
    
    Ey = cos(gama1*a + phi)*exp( gama2*(x+a)) .* (x < -a)+...
        cos(gama1*x - phi)                    .* (x >=-a & x <= a)+...
        cos(gama1*a - phi)*exp(-gama2*(x-a))  .* (x > a);
    
    subplot(2,2,count);
    plot(x*1e6,Ey);
    xlim([-3e6*a 3e6*a]);
    ylim([-1 1]);
    count = count + 1;
end
figure(2);plot(neff, 'ro--');