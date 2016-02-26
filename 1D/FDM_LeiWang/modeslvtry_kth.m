I = sqrt(-1);
Pi = 3.1416;
lambda = 1.55e-6; %1.15;
k0 = 2*Pi/lambda;


dx = 0.1e-9;
CldWid = 1.85e-6;
SCHWid1 = 0.025e-6;
SCHWid2 = 0.023e-6;
QwWid = 8.8e-9;
BrWid = 9.7e-9;
StpWid = 10e-9;
CoerWid = 0.4e-6;
SubWid = 4e-6;


CldInd = 3.1694;
CldAct = 3.1694;
SCHInd1 = 3.3766;
SCHInd2 = 3.2891;
QwInd = 3.59;
BrInd = 3.3766;
StpInd = 3.1694;
CoerInd = 3.37;
SubInd = 3.1694;

vThickness = [SubWid,CoerWid,StpWid,SCHWid2,SCHWid1,BrWid,QwWid,BrWid,QwWid,BrWid,QwWid,BrWid,QwWid,BrWid,SCHWid1,SCHWid2,CldWid];
vNeffAct =   [SubInd,CoerInd,StpInd,SCHInd2,SCHInd1,BrInd,QwInd,BrInd,QwInd,BrInd,QwInd,BrInd,QwInd,BrInd,SCHInd1,SCHInd2,CldAct];
vIsAct  =     [0,    0,        0,     0,     0,      1,    1,     1,   1,    1,     1,    1,    1,    1,    0,    0,     0];

vNeffPas =   [SubInd,CoerInd,StpInd,CldInd,CldInd,CldInd,CldInd,CldInd,CldInd,CldInd,CldInd,CldInd,CldInd,CldInd,CldInd,CldInd,CldInd];
n = length(vNeffAct);
nLayers = zeros(1,n);
iLayers = ones(1,n+1);
eepsAct = [];
eepsPas = [];
vIsAct1 = [];
for i=1:n
    nLayers(i) = round(vThickness(i)/dx);
    iLayers(i+1) = iLayers(i) + nLayers(i);
    eepsTemp = vNeffAct(i)^2*ones(1,nLayers(i));
    eepsAct = [eepsAct,eepsTemp];
    vIsAct1 = [vIsAct1,vIsAct(i)*ones(1,nLayers(i))];
    
    eepsTemp = vNeffPas(i)^2*ones(1,nLayers(i));
    eepsPas = [eepsPas,eepsTemp];
end;

% viS = zeros(1,7);
% viE = zeros(1,7);
% for i = 1:7
%     viS(i) = iLayers(2*i+2);
%     viE(i) = iLayers(2*i+3)-1;
% end;

nx = length(eepsAct)-1;
vIsAct1 = vIsAct1(1:nx);
%boundary condition is zero boundary
boundary = 'oo';  

%the number of modes to be caculated
nmodes = 1;
% the initial guess of the mode index
guess = 4;
% caculate the mode
[vAct dAct] = Modesolver1D(eepsAct,k0,dx,nmodes,guess,lambda,boundary,'TE');
[vPas dPas] = Modesolver1D(eepsPas,k0,dx,nmodes,guess,lambda,boundary,'TE');

neff = zeros(1,nmodes);
Gama = zeros(1,nmodes);
Intg = zeros(1,nmodes);
for(i=1:nmodes)
    neff(i) = sqrt(dAct(i,i))/k0
    
    Gama(i) = sum(abs(vAct(:,i).^2).*vIsAct1.')/sum(abs(vAct(:,i).^2));
    vx = dx*(1:nx)*1e6;
    Intg(i) = Intergral(vx,vAct(:,i),vx,vPas(:,i),dx)
     figure;
     hold;
      plot(vx,( abs(vAct(:,i)))/max( abs(vAct(:,i))),'r');
      plot(vx,( abs(vPas(:,i)))/max( abs(vPas(:,i))),'b');
     hold;
%      figure;
%       plot(dx*(1:nx)*1e6,angle(vAct(:,i))*180/Pi);
end;

%д���ļ�
% v(:,1) = v(:,1)/sqrt(max(abs(v(:,1)).^2));
% fid=fopen('d:\rsoft\examples\tmspmode1.m00','w');
%     
%         fprintf(fid,'/rn,a,b/nx0/ls1 \r\n');
%         fprintf(fid,'%d %g %g %g OUTPUT_REAL_IMAG %g %g \r\n', nx, -(a/2+tair)*1e6,((a/2)+twide-tair-a)*1e6,0, real(neff(1)),imag(neff(1)));
% 
%  fprintf(fid,'%g %g \r\n',[real(v(:,1))';imag(v(:,1))']);
%  fclose(fid);
 
%  fid=fopen('e:\resultsp\result1.mon','r');
%     for(j=1:5)
%         tline = fgetl(fid);
%     end;
%     temp = fscanf(fid,'%g %g %g',[3 inf]);
%     hold on;
%     plot(temp(1,:),(10*log10(temp(2,:))),'k');