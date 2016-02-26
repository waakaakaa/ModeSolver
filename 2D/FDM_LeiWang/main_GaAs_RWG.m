clear
clc
close all

I = sqrt(-1);
Pi = 3.1416;
lambda = 0.85e-6; %1.15;
k0 = 2*Pi/lambda;

thnCld = 1.0e-6;
thnQw = 0.107e-6;
thnSub = 2e-6;

CldInd = 3.09;
QwInd = 3.4427;

SubInd = 3.09;



vThickness = [thnSub,thnQw,0.1e-6,thnCld,0.5e-6];
vNeffAct =   [SubInd,QwInd, CldInd,CldInd,1.6];
vIsAct  =     [0,        0,     0, 1,  0, ];
iStart= 2;
iEnd = 2;

twide = 9e-6;
tthn = 5e-6;
ridW = 1e-6;

dx = 0.02e-6;
dy = 0.01e-6;
nmodes = 1;

vw = [1,2,3,9]*1e-6;


n = length(vw);

vGama0 = zeros(1,n);%基模限制因子
vGama1 = zeros(1,n);%一阶模限制因子
vneff = zeros(1,n);
for i=1:n
    i
    %------------------------获得折射率分布矩阵---------------%
    [eeps,nx,ny,mStEnd] = GenerateEps(vThickness,vNeffAct,[twide,twide,twide,vw(i),twide],1.55,twide,dx,dy);%取得左右对称的多层波导结构的折射率分布矩阵

    guess = QwInd;
% caculate the mode
    boundary = 'ooss';  % the boundary condition must be set 'ss' for slab
    field = 'ex';
    A = svbuildmtx (lambda,dx,dy,eeps,boundary,field);% the difference 
    
    bIsDrawing = 1;
    [mE0,neff,vcellE] = sveigenmodes(A,guess,nmodes,lambda,nx,ny,[],[],bIsDrawing);
    vneff(i) = neff(1);
    
    vGama0(i) = GetConfineBylayerNumbers(vcellE{1},nx,ny,iStart,iEnd,mStEnd);%获得iStart layer 到iEnd layer之间光能量占总能量的比例
    %vGama1(i) = GetConfineBylayerNumbers(vcellE{2},nx,ny,iStart,iEnd,mStEnd);%获得iStart layer 到iEnd layer之间光能量占总能量的比例
end
%close all;
figure
plot(vw*1e6,(vneff-min(vneff))/min(vneff)*lambda*1e9);
figure;
hold on;
plot(vw*1e6,vGama0,'r');
%plot(vw*1e6,vGama1,'g');
% [filename, pathname, filterindex] = uiputfile( ...
% {'*.mat';},'Save as');
% str = [pathname,filename];
% save (str)





