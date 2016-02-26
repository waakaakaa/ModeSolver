clc;
close all;
clear;
I = sqrt(-1);
Pi = 3.1416;
lambda0 = 0.45e-6; 
k0 = 2*Pi/lambda0;

dx = 0.5e-9;

pConWid = 0.1e-6;
nConWid = 4.5e-6;
pSLSWid = 0.6e-6;
nSLSWid = 1e-6;
pWgWid = 0.07e-6;
nWgWid = 0.07e-6;
EBLWid = 0.015e-6;
LBWid = 0.02e-6;
QwWid = 0.003e-6;
BrWid = 0.015e-6;
SubWid = 1e-6;

pConInd = 2.469;
nConInd = 2.469;
pSLSInd = 2.4305;
nSLSInd = 2.4305;
pWgInd = 2.531;
nWgInd = 2.531;
EBLInd = 2.418;
LBInd = 2.493;
QwInd = 2.733;
BrInd = 2.493;
SubInd = 1.783;

vThickness = [SubWid,nConWid,nSLSWid,nWgWid,BrWid,QwWid,BrWid,QwWid,BrWid,QwWid,LBWid,EBLWid,pWgWid,pSLSWid,pConWid,1e-6];
vNeffAct = [SubInd,nConInd,nSLSInd,nWgInd,BrInd,QwInd,BrInd,QwInd,BrInd,QwInd,LBInd,EBLInd,pWgInd,pSLSInd,pConInd,1];

[eeps] = GetEpslProf1D(vThickness,vNeffAct,0,dx);%epsl profile in waveguide section 

boundary = 'oo';  
    %the number of modes to be caculated
    nmodes = 1;
    % the initial guess of the mode index
    guess = 2.6;
[vEte d] = Modesolver1D(eeps,k0,dx,nmodes,guess,lambda0,boundary,'TE');


n = length(vThickness);
vn = round(vThickness/dx);
vInd = [];
 for i = 1:n;
        vInd = [vInd,vNeffAct(i)*ones(1,vn(i))];
 end
 
figure;
[ax,h1,h2] = plotyy(dx*(1:sum(vn)-1)*1e6,abs(vEte.^2),dx*(1:sum(vn)-1)*1e6,vInd(1:sum(vn)-1));
xlabel(ax(1),'distance from substrate(\mum)');ylabel(ax(1),'Intensity(a.u.)')       
xlabel(ax(2),'distance from substrate(\mum)');ylabel(ax(2),'Index') 
hold on
neff = sqrt(d)/k0


