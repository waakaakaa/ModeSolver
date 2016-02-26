function [ Ex ] = plotExMultiLayer( neff ,x)

global lambda0;
global n0;
global refractiveIndex;
global nLast;
global thickness;

n = length(refractiveIndex);
k0 = 2*pi/lambda0;

gama0 = k0 * sqrt(neff.^2 - n0^2);
gamaLast = k0 * sqrt(neff.^2 - nLast^2);
gama = k0 * sqrt(abs(neff.^2 - refractiveIndex.^2));


T = zeros(1,n);
A = zeros(1,n);
B = zeros(1,n);

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
ALast = (A(end)*cos(gama(end)*thickness(end)) + B(end)*sin(gama(end)*thickness(end))) .* (neff<=refractiveIndex(end))...
    +(A(end)*cosh(gama(end)*thickness(end)) + B(end)*sinh(gama(end)*thickness(end))) .* (neff>refractiveIndex(end));

sigmaB = zeros(1,n + 1);
temp = 0;
sigmaB(1) = 0;
for t = 1:n
    temp = temp + thickness(t);
    sigmaB(t+1) = temp;
end

Ex = exp(gama0*x).*(x<0);
for t = 1:n
    ExTemp = ( A(t)* cos(gama(t)*(x - sigmaB(t))) + B(t)* sin(gama(t)*(x - sigmaB(t))) ).*(neff< refractiveIndex(t))...
        + ( A(t)*cosh(gama(t)*(x - sigmaB(t))) + B(t)*sinh(gama(t)*(x - sigmaB(t))) ).*(neff>=refractiveIndex(t));
    ExTemp = ExTemp .* (x>=sigmaB(t) & x<sigmaB(t+1));
    Ex = Ex + ExTemp;
end
Ex = Ex + ALast*exp(-gamaLast*(x - sigmaB(end))).*(x>=sigmaB(end));

end