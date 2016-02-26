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
thnBuf = 0.05e-6;
thnGrin  = 0.01e-6;


CldInd = 3.1694;
CldAct = 3.1694;
QwInd = 3.43;
GrinInd = 3.3;
CoreInd = 3.37;
SubInd = 3.1694;


vThickness = [thnSub,0.12e-6,thnQw,thnGrin, thnBuf, thnCore,thnCld,0.2e-6];
vNeffAct =   [SubInd,GrinInd,QwInd,GrinInd,CldAct, CoreInd, CldAct,1.55];
vIsAct  =     [0,        0,     1,  0,         0,    0,   0,    0];
iStart= 3;
iEnd = 3;

twide = 8e-6;
tthn = 5e-6;
ridW = 4e-6;

dx = 0.02e-6;
dy = 0.01e-6;
nmodes = 3;

vw = [0.8:0.2:1.6]*1e-6;


n = length(vw);

vGama0 = zeros(1,n);%��ģ��������
vGama1 = zeros(1,n);%һ��ģ��������
vneff = zeros(1,n);
for i=1:n
    i
    %------------------------��������ʷֲ�����---------------%
    [eeps,nx,ny,mStEnd] = GenerateEps(vThickness,vNeffAct,[twide,twide,twide,twide,twide,vw(i),2e-6+vw(i)-vw(1),twide],1.55,twide,dx,dy);%ȡ�����ҶԳƵĶ�㲨���ṹ�������ʷֲ�����

    guess = CoreInd;
% caculate the mode
    boundary = 'oooo';  % the boundary condition must be set 'ss' for slab
    field = 'ex';
    A = svbuildmtx (lambda,dx,dy,eeps,boundary,field);% the difference 
    
    bIsDrawing = 1;
    [mE0,neff,vcellE] = sveigenmodes(A,guess,nmodes,lambda,nx,ny,[],[],bIsDrawing);
    vneff(i) = neff(1);
    
    vGama0(i) = GetConfineBylayerNumbers(vcellE{1},nx,ny,iStart,iEnd,mStEnd);%���iStart layer ��iEnd layer֮�������ռ�������ı���
    vGama1(i) = GetConfineBylayerNumbers(vcellE{2},nx,ny,iStart,iEnd,mStEnd);%���iStart layer ��iEnd layer֮�������ռ�������ı���
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





