clear
clc
close all

I = sqrt(-1);
Pi = 3.1416;
lambda = 1.55e-6; %1.15;
k0 = 2*Pi/lambda;

thnCld = 1.8e-6;
thnCld2 = 0.15e-6;
thnCore = 0.4e-6;
thnSub = 2e-6;
thnEtc = 0.01e-6;


CldInd = 3.1694;

CoreInd = 3.43;
SubInd = 3.1694;


vThickness = [thnSub,thnCore,thnCld2,thnCld,0.2e-6,1e-6];
vNeff =   [SubInd,CoreInd,CldInd,CldInd,CoreInd ,1.0];


twide = 16e-6;%���������ܿ���
tthn = 5e-6;
ridW = 2e-6;%��������

dx = 0.02e-6;
dy = 0.02e-6;
nmodes = 3;

vw = [1:0.5:3]*1e-6;%��������ɨ������


n = length(vw);

vneff0 = zeros(1,n);%��ģ������
vneff1 = zeros(1,n);%1��ģ������
vneff2 = zeros(1,n);%2��ģ������
for i=1:n
    i
    %------------------------��������ʷֲ�����---------------%
    [eeps,nx,ny,mStEnd] = GenerateEps(vThickness,vNeff,[twide,twide,twide,vw(i),vw(i),twide],1.0,twide,dx,dy);%ȡ�����ҶԳƵĶ�㲨���ṹ�������ʷֲ�����

    guess = CoreInd;
% caculate the mode
    boundary = 'oooo';  % the boundary condition must be set 'ss' for slab
    field = 'ey';
    A = svbuildmtx (lambda,dx,dy,eeps,boundary,field);% the difference 
    
    bIsDrawing = 1;
    [mE0,neff,vcellE] = sveigenmodes(A,guess,nmodes,lambda,nx,ny,[],[],bIsDrawing);
    vneff0(i) = neff(1);
    vneff1(i) = neff(2);
    vneff2(i) = neff(3);
end
%close all;
figure;
plot(vw*1e6,vneff0,'r.');
hold
plot(vw*1e6,vneff1,'g.');
plot(vw*1e6,vneff2,'b.');

% [filename, pathname, filterindex] = uiputfile( ...
% {'*.mat';},'Save as');
% str = [pathname,filename];
% save (str)




