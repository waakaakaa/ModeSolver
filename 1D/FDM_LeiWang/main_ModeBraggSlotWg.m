I = sqrt(-1);
Pi = 3.1416;
lambda0 = 1.55e-6; %1.15;
k0 = 2*Pi/lambda0;


dx = 0.5e-9;
nSiO2 = 1.44;
nSi = 3.445;

thnSiO2 = 1/4*lambda0/nSiO2;
thnSi = 1/32*lambda0/nSi;
np = 8;

vNeff = zeros(1,np*2+2);
vThickness = zeros(1,np*2+2);
for i=1:np
    vNeff(i*2) = nSi;
    vThickness(i*2) = thnSi;
    vNeff(i*2+1) = nSiO2;
    vThickness(i*2+1) = thnSiO2;
end
%vL((nn-1)/2+1) = d1*2;
vThickness(np+1) = thnSiO2*2;
vNeff(np+1) = 1.55;
vThickness(1) = 2e-6;
vNeff(1) = nSiO2/1;
vThickness(np*2+2) = 2e-6;
vNeff(np*2+2) = nSiO2/1;

n = length(vNeff);
nLayers = zeros(1,n);
iLayers = ones(1,n+1);
eeps = [];

for i=1:n
    nLayers(i) = round(vThickness(i)/dx);
    iLayers(i+1) = iLayers(i) + nLayers(i);
    eepsTemp = vNeff(i)^2*ones(1,nLayers(i));
    eeps = [eeps,eepsTemp];
end

nx = length(eeps)-1;
vx = [0:nx-1]*dx;

boundary = 'oo';  
%the number of modes to be caculated
nmodes = 1;
% the initial guess of the mode index
guess = 4;
% caculate the mode
[vE d] = Modesolver1D(eeps,k0,dx,nmodes,guess,lambda0,boundary,'TE');
figure;
plot(vx,abs(vE.^2)/max(abs(vE.^2)));
hold on;
plot(vx,sqrt(eeps(1:nx))/max(sqrt(eeps(1:nx))),'r');
