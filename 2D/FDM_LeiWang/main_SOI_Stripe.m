clear
clc
close all
I = sqrt(-1);
Pi = 3.1416;

nSub = 1.44;
nCore = 3.48;
nCld = 1.331;

thnSub = 2.0e-6;
thnCore = 0.22e-6;
thnCld = 2e-6;

twide = 8e-6;
tthn = 5e-6;
ridW = 0.5e-6;

dx = 0.02e-6;
dy = 0.01e-6;
nmodes = 1;
%------------------------获得折射率分布矩阵---------------%
[eeps,nx,ny,mStEnd] = GenerateEps([thnCld,thnCore,thnSub],[nCld,nCore,nSub],[twide,ridW,twide],nCld,twide,dx,dy);%取得左右对称的多层波导结构的折射率分布矩阵

vlambda = ([-30:30]+1550)*1e-9;
n = length(vlambda);
vneff = zeros(1,n);
for i=1:n
    i
    lambda = vlambda(i);
    guess = nCore;
% caculate the mode
    boundary = 'oooo';  % the boundary condition must be set 'ss' for slab
    field = 'ey';
    A = svbuildmtx (lambda,dx,dy,eeps,boundary,field);% the difference 
    
    bIsDrawing = 0;
    [g,neff] = sveigenmodes(A,guess,nmodes,lambda,nx,ny,[],[],bIsDrawing);
    vneff(i) = neff;
end
[p] = polyfit(vlambda,vneff,1);
a0 = p(2);
a1 = p(1);
[filename, pathname, filterindex] = uiputfile( ...
{'*.mat';},'Save as');
str = [pathname,filename];
save (str)





