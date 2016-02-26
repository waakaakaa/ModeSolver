function [ N ] = plotNMultiLayer(x)

global n0;
global refractiveIndex;
global nLast;
global thickness;

n = length(refractiveIndex);

sigmaB = zeros(1,n + 1);
temp = 0;
sigmaB(1) = 0;
for t = 1:n
    temp = temp + thickness(t);
    sigmaB(t+1) = temp;
end

N = n0.*(x<0);
for t = 1:n
    Ntemp = refractiveIndex(t);
    Ntemp = Ntemp .* (x>=sigmaB(t) & x<sigmaB(t+1));
    N = N + Ntemp;
end
N = N + nLast.*(x>=sigmaB(end));

end