clc;
clear;

global lambda0;
lambda0 = 1.55e-6;
k0 = 2*pi/lambda0;

global n1;
n1 = 3.7262;

global n2;
n2 = 3.5214;

global n3;
n3 = 3.5941;

global a;
a = 0.5e-6;

global d;
d = 0.5e-6;

global m;
m = 0;



dd = (0: 0.1: 2.0)*1e-6;
aa = (0.0: 0.05: 0.4)*1e-6;
neff = zeros(1, length(dd));
count = 1;
figure(1);
hold on;
box on;
for a = aa
    for d = dd
        %         options = optimoptions('fsolve','Display','iter');
        options = optimset('Display','none');
        [x,fval] = fsolve(@calBeta5LayerSymmetric, k0*(n1-0.01), options);
        
        beta0 = real(x(1));
        neff(1, count) = beta0/k0;
        count = count + 1;
    end
    plot(dd*1e6, neff(1,:), 'r.-');
    count = 1;
end
xlim([0 2]);
ylim([3.595 3.685]);
hold off;



dd = [0 0.2 1.0]*1e-6;
aa = (0: 0.1: 1.5)*1e-6;
neff = zeros(1, length(aa));
count = 1;
figure(2);
hold on;
box on;
for m = [0 1 2 3]
    for d = dd
        for a = aa
            %             options = optimoptions('fsolve','Display','iter');
            options = optimset('Display','none');
            [x,fval] = fsolve(@calBeta5LayerSymmetric, k0*(n2+0.01), options);
            
            beta0 = real(x(1));
            neff(1, count) = beta0/k0;
            count = count + 1;
        end
        plot(aa*1e6, neff(1,:), 'r-');
        count = 1;
    end
end
xlim([0 1.5]);
ylim([3.595 3.73]);
hold off;