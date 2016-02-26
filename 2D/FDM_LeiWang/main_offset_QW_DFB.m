clear
clc
close all

I = sqrt(-1);
Pi = 3.1416;
lambda = 1.55e-6; %1.15;
k0 = 2*Pi/lambda;

thnCld = 1.8e-6;
thnQw = 0.10e-6;
thnCore = 0.4e-6;
thnSub = 3e-6;
thnEtc = 0.01e-6;
thnBuf = 0.01e-6;
thnGrin  = 0.01e-6;


CldInd = 3.1694;
CldAct = 3.1694;
QwInd = 3.43;
GrinInd = 3.3;
CoreInd = 3.37;
SubInd = 3.1694;


vThickness = [thnSub,thnCore,thnBuf,thnQw,0.05e-6,  thnCld,0.2e-6];
vNeffAct =   [SubInd,CoreInd,CldAct,QwInd,GrinInd,  CldAct,1.55];
vIsAct  =     [0,        0,     1,  0,         0,       0,    0];
iStart= 4;
iEnd = 4;

twide = 8e-6;
tthn = 5e-6;
ridW = 4e-6;

dx = 0.02e-6;
dy = 0.01e-6;
nmodes = 3;

vw = [0.8:0.2:1.6]*1e-6;


n = length(vw);

vGama0 = zeros(1,n);%基模限制因子
vGama1 = zeros(1,n);%一阶模限制因子
vneff = zeros(1,n);
for i=1:n
    i
    %------------------------获得折射率分布矩阵---------------%
    [eeps,nx,ny,mStEnd] = GenerateEps(vThickness,vNeffAct,[twide,vw(i),twide,twide,twide,twide,twide,twide],1.55,twide,dx,dy);%取得左右对称的多层波导结构的折射率分布矩阵

    guess = CoreInd;
% caculate the mode
    boundary = 'oooo';  % the boundary condition must be set 'ss' for slab
    field = 'ex';
    A = svbuildmtx (lambda,dx,dy,eeps,boundary,field);% the difference 
    
    bIsDrawing = 1;
    [mE0,neff,vcellE] = sveigenmodes(A,guess,nmodes,lambda,nx,ny,[],[],bIsDrawing);
    vneff(i) = neff(1);
    
    vGama0(i) = GetConfineBylayerNumbers(vcellE{1},nx,ny,iStart,iEnd,mStEnd);%获得iStart layer 到iEnd layer之间光能量占总能量的比例
    vGama1(i) = GetConfineBylayerNumbers(vcellE{2},nx,ny,iStart,iEnd,mStEnd);%获得iStart layer 到iEnd layer之间光能量占总能量的比例
end
%close all;
figure
plot(vw*1e6,(vneff-min(vneff))/min(vneff)*lambda*1e9);
figure;
hold on;
plot(vw*1e6,vGama0,'r');
plot(vw*1e6,vGama1,'g');
% [filename, pathname, filterindex] = uiputfile( ...
% {'*.mat';},'Save as');
% str = [pathname,filename];
% save (str)





