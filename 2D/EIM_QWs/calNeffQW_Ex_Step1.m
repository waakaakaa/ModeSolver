function [ F ] = calNeffQW_Ex_Step1( neff )

global lambda0;
global n0;
global refractiveIndex;
global nLast;
global thickness;
global nn;

n = length(refractiveIndex);
k0 = 2*pi/lambda0;

gama0 = k0 * sqrt(neff.^2 - n0^2);
gamaLast = k0 * sqrt(neff.^2 - nLast^2);


T = zeros(1,n);
A = zeros(1,n);
B = zeros(1,n);
gama = k0 * sqrt(abs(neff.^2 - refractiveIndex.^2));


for t = (1:n)
    if t == n
        T(t) = gama(t)/gamaLast;
    else
        T(t) = gama(t)/gama(t+1);
    end
    if t == 1
        A(t) = 1;
        B(t) = gama0 / gama(1);
    else
        A(t) = (A(t-1)*cos(gama(t-1)*thickness(t-1)) + B(t-1)*sin(gama(t-1)*thickness(t-1))) .* (neff<=refractiveIndex(t-1))...
            +(A(t-1)*cosh(gama(t-1)*thickness(t-1)) + B(t-1)*sinh(gama(t-1)*thickness(t-1))) .* (neff>refractiveIndex(t-1));
        B(t) = T(t-1)*(-A(t-1)*sin(gama(t-1)*thickness(t-1)) + B(t-1)*cos(gama(t-1)*thickness(t-1))) .* (neff<=refractiveIndex(t-1))...
            +T(t-1)*(A(t-1)*sinh(gama(t-1)*thickness(t-1)) + B(t-1)*cosh(gama(t-1)*thickness(t-1))) .* (neff>refractiveIndex(t-1));
    end
end

% This equation is only correct for 0th mode.
F = (    gama(end)*thickness(end) - nn*pi - atan( (A(end) + B(end)*T(end))/(A(end)*T(end) - B(end)) )    ) .* (neff<refractiveIndex(end))...
    + (  tanh(gama(end)*thickness(end))  +  (A(end) + B(end)*T(end))/(A(end)*T(end) + B(end)) )     .* (neff>=refractiveIndex(end));


end