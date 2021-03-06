clc
clear
I = sqrt(-1);
Pi = 3.1416;

dx=0.02e-6;
dy=0.001e-6;
totalwid = 12e-6;
mesawid = 10e-6;
rgwid = 3e-6;
%waveguidewid=20e-6;

ncld=1;

pConTh = 0.1e-6;
nConTh1 = 2e-6;         % n_GaN设为两层，具有不同的宽度
nConTh2 = 1e-6;
pSLSTh = 0.6e-6;
nSLSTh = 1e-6;
pWgTh = 0.06e-6;
nWgTh = 0.06e-6;
EBLTh = 0.015e-6;
LBTh = 0.02e-6;
QwTh = 0.003e-6;
BrTh = 0.01e-6;
SubTh = 1e-6;

pConInd = 2.469;
nConInd1 = 2.469;
nConInd2 = 2.469;
pSLSInd = 2.42;
nSLSInd = 2.42;
pWgInd = 2.506;
nWgInd = 2.506;
EBLInd = 2.377;
LBInd = 2.506;
QwInd = 2.733;
BrInd = 2.506;
SubInd = 1.783;

pConWid = rgwid;
nConWid1 = totalwid;
nConWid2 = mesawid;
pSLSWid = rgwid;
nSLSWid = mesawid;
pWgWid = mesawid;
nWgWid = mesawid;
EBLWid = mesawid;
LBWid = mesawid;
QwWid = mesawid;
BrWid = mesawid;
SubWid = totalwid;

vThickness = [SubTh,nConTh1,nConTh2,nSLSTh,nWgTh,BrTh,QwTh,BrTh,QwTh,BrTh,QwTh,BrTh,QwTh,LBTh,EBLTh,pWgTh,pSLSTh,pConTh,1e-6];
vNeffact = [SubInd,nConInd1,nConInd2,nSLSInd,nWgInd,BrInd,QwInd,BrInd,QwInd,BrInd,QwInd,BrInd,QwInd,LBInd,EBLInd,pWgInd,pSLSInd,pConInd,1];
vWidth = [SubWid,nConWid1,nConWid2,nSLSWid,nWgWid,BrWid,QwWid,BrWid,QwWid,BrWid,QwWid,BrWid,QwWid,LBWid,EBLWid,pWgWid,pSLSWid,pConWid,totalwid];


%-----------% 八量子阱片%
% vThickness=[1.5,0.025,0.025,0.025,0.025,tqw,tbarrier,tqw,tbarrier,tqw,tbarrier,tqw,tbarrier,tqw,tbarrier,tqw,tbarrier,tqw,tbarrier,tqw,0.025,0.025,0.025,0.15,1.5]*1e-6;
% vNeffact=[3.1694,3.2525,3.2829,3.31,3.364,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.364,3.31,3.2525,3.1694,3.1694];
% kk=length(vNeffact);
% vWidth=[];
% for i=1:(kk-2)
%     vWidth=[vWidth,totalwid];
% end
% vWidth=[vWidth,3e-6,3e-6];

[eeps,nx,ny,mStEnd] = GenerateEps(vThickness,vNeffact,vWidth,ncld,totalwid,dx,dy);

lambda = 0.45e-6;
    guess = 2.5;
    % caculate the mode
    boundary = 'oooo';  % 参svbuiltmx里面的说明
    field = 'ey';
    A = svbuildmtx (lambda,dx,dy,eeps,boundary,field);
     bIsDrawing = 1;
     nmodes=1;
    [g,neff] = sveigenmodes(A,guess,nmodes,lambda,nx,ny,[],[],bIsDrawing);
 neff


