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

n1 = 1.5;
n2 = 1;
n3 = 1;
n4 = 1;
n5 = 1;
lambda0 = 1.55e-6;
n = 0;

k0 = 2*pi/lambda0;

figure(1);hold on;box on;
B_vector = linspace(0, 4, 50);
count = 1;
for m = [0 1 2]
    P = zeros(1,length(B_vector));
    for B = B_vector
        b = B*lambda0/2./sqrt(n1^2-n2^2);
        a = b*2;
        %         options = optimoptions('fsolve','Display','final');
        options = optimset('Display','none');
        
        [kx,fval1] = fsolve(@calKx_Ey, 5e5, options);
        [ky,fval2] = fsolve(@calKy_Ey, 5e5, options);
        
        kz = sqrt(k0^2*n1^2 - kx^2 - ky^2);
        neff = kz/k0;
        if neff > n2
            P(1,count) = (neff^2 - n2^2)/(n1^2 - n2^2);
        end
        count = count + 1;
    end
    plot(B_vector,P,'b.-')
    count = 1;
end

for m = [0 1 2]
    P = zeros(1,length(B_vector));
    for B = B_vector
        b = B*lambda0/2./sqrt(n1^2-n2^2);
        a = b*2;
        neff = sqrt(...
            n1^2 - ...
            (m+1)^2*pi^2 / ...
            (k0*a + sqrt(n1^2 - n3^2) + sqrt(n1^2 - n5^2))^2 -...
            (n+1)^2*pi^2 / ...
            (k0*b + (n2^2/n1^2)*sqrt(n1^2 - n2^2) + (n4^2/n1^2)*sqrt(n1^2 - n4^2))^2  );
        if neff > n2
            P(1,count) = (neff^2 - n2^2)/(n1^2 - n2^2);
        end
        count = count + 1;
    end
    plot(B_vector,P,'r.-')
    count = 1;
end
hold off;