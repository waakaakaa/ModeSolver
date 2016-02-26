clear
clc
close all
I = sqrt(-1);
Pi = 3.1416;

nSub = 3.1694;
nCore = 3.37;
nCld = 3.1694;
nCov = 1.46;

thnSub = 2.0e-6;
thnCore = 0.4e-6;
thnCld = 1e-6;

twide = 18e-6;
tthn = 5e-6;
ridW = 3e-6;

dx = 0.02e-6;
dy = 0.02e-6;
nmodes = 4;
%------------------------获得折射率分布矩阵---------------%
[eeps,nx,ny,mStEnd] = GenerateEps([0.5e-6,thnCld-0.15e-6,0.15e-6,thnCore,thnSub],[nCov,nCld,nCld,nCore,nSub],[twide,ridW,twide,twide,twide,twide],nCov,twide,dx,dy);%取得左右对称的多层波导结构的折射率分布矩阵

lambda = 1.55e-6;
    guess = nCore;
% caculate the mode
    boundary = 'oooo';  % the boundary condition must be set 'ss' for slab
    field = 'ey';
    A = svbuildmtx (lambda,dx,dy,eeps,boundary,field);% the difference 
    
    bIsDrawing = 1;
    [g,neff] = sveigenmodes(A,guess,nmodes,lambda,nx,ny,[],[],bIsDrawing);
    





