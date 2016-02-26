clc;
close all;
clear;
I = sqrt(-1);
Pi = 3.1416;
lambda0 = 1.55e-6; %1.15;
k0 = 2*Pi/lambda0;


dx = 0.5e-9;
nCld = 3.1694;
nCore = 3.2219;

thnCld = 1.0e-6;
thnCore = 0.6e-6;
nAir = 1.44;


vEtdep = [0,[0.5:0.025:1.0]*1e-6];%刻蚀深度
nn = length(vEtdep);
vneffte = zeros(1,nn);
vnefftm = zeros(1,nn);
vCte = zeros(1,nn);%耦合系数
vCtm = zeros(1,nn);%耦合系数

for i=1:nn
    vEtdep(i)
    vThn = [1e-6,vEtdep(i),thnCld-vEtdep(i),thnCore,6e-6];
    vN1 = [nAir,   nAir        nCld,       nCore,nCld];
    [eeps] = GetEpslProf1D(vThn,vN1,0,dx);%epsl profile in waveguide section 
    vx = [0:length(eeps)-2]*dx;

    boundary = 'oo';  
    %the number of modes to be caculated
    nmodes = 1;
    % the initial guess of the mode index
    guess = 4;
    % caculate the mode
    [vEtm d] = Modesolver1D(eeps,k0,dx,nmodes,guess,lambda0,boundary,'TM');
    vnefftm(i) = sqrt(d)/k0;
    [vEte d] = Modesolver1D(eeps,k0,dx,nmodes,guess,lambda0,boundary,'TE');
    vneffte(i) = sqrt(d)/k0;
    if(i==1)
        vEte0 = vEte;
        vEtm0 = vEtm;
    end
    vCte(i) = abs(IntergralwithPhase(vx,vEte,vx,vEte0,dx)).^2/sum(abs(vEte.^2))/sum(abs(vEte0.^2));
    vCtm(i) = abs(IntergralwithPhase(vx,vEtm,vx,vEtm0,dx)).^2/sum(abs(vEtm.^2))/sum(abs(vEtm0.^2));
end
vgama = (vneffte(1)-vnefftm(1))./((vneffte(2:nn)-vnefftm(2:nn))-(vneffte(1)-vnefftm(1)));
figure;
plot(vEtdep(2:nn),vneffte(2:nn)-vnefftm(2:nn),'.');
figure;
subplot(3,1,1);
plot(vEtdep(2:nn),vgama,'-');
subplot(3,1,2);
plot(vEtdep(2:nn),vCte(2:nn),'b');
hold on;
plot(vEtdep(2:nn),vCtm(2:nn),'r');
subplot(3,1,3);
plot(vgama,vCte(2:nn),'b');
hold on;
plot(vgama,vCtm(2:nn),'r');

% figure;
% plot(vx,abs(vE.^2)/max(abs(vE.^2)));

