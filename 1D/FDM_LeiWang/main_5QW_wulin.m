clc;
close all;
clear;
I = sqrt(-1);
Pi = 3.1416;
lambda0 = 1.55e-6; %1.15;
k0 = 2*Pi/lambda0;


dx = 0.5e-9;
Cld1Wid = 0.15e-6;
Cld2Wid = 1.5e-6;
SCHWid = 0.025e-6;

QwWid = 5.5e-9;
BrWid = 10e-9;
StpWid = 10e-9;
SubWid = 10e-6;


CldInd = 3.1694;
CldAct = 3.1694;
SCHInd = 3.37;
QwInd = 3.59;
BrInd = 3.37;
StpInd = 3.1694;
CoerInd = 3.37;
SubInd = 3.1694;

vThickness = [SubWid,SCHWid,SCHWid,SCHWid,SCHWid,QwWid,BrWid,QwWid,BrWid,QwWid,BrWid,QwWid,BrWid,QwWid,SCHWid,SCHWid,SCHWid,Cld1Wid,Cld2Wid,1e-6];
vNeffAct =   [SubInd,3.2525,3.2829,3.31,  3.364, QwInd,BrInd,QwInd,BrInd,QwInd,BrInd,QwInd,BrInd,QwInd,3.364, 3.31,  3.2525,CldAct,CldAct,1.55];
vNeffAct1 =   [SubInd,3.2525,3.2829,3.31,  3.364, QwInd,BrInd,QwInd,BrInd,QwInd,BrInd,QwInd,BrInd,QwInd,3.364, 3.31,  3.2525,CldAct,1.55,1.55];
[eeps] = GetEpslProf1D(vThickness,vNeffAct,0,dx);%epsl profile in waveguide section 

boundary = 'oo';  
    %the number of modes to be caculated
    nmodes = 1;
    % the initial guess of the mode index
    guess = 4;
[vEte d] = Modesolver1D(eeps,k0,dx,nmodes,guess,lambda0,boundary,'TE');
figure;
plot(abs(vEte.^2));
neff = sqrt(d)/k0


