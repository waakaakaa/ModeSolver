clc
clear
I = sqrt(-1);
Pi = 3.1416;

dx=0.02e-6;
dy=0.002e-6;
totalwid=12e-6;
%waveguidewid=20e-6;

ncld=1;

% vThickness=[1.5,0.15,0.02,0.02,0.01,0.01,0.0055,0.01,0.0055,0.01,0.0055,0.01,0.0055,0.01,0.0055,0.01,0.01,0.02,0.02,5]*1e-6;
% vNeffact=[1.55,3.1694,3.2525,3.31,3.364,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.364,3.31,3.2829,3.1694];
% vWidth=[totalwid,totalwid];
%
% ll=length(vThickness);
% for k=1:1:ll
%     temp(k)=vThickness(ll+1-k);
% end
% vThickness=temp;
%
% for k=1:1:ll
%     temp(k)=vNeffact(ll+1-k);
% end
% vNeffact=temp;
%
% for i=1:1:(ll-2)
%     vWidth=[vWidth,totalwid];
% end
%background

tbarrier=0.01;
tqw=0.0055;

vThickness=[1.5,0.02,0.02,0.01,tbarrier,tqw,tbarrier,tqw,tbarrier,tqw,tbarrier,tqw,tbarrier,tqw,tbarrier,0.01,0.02,0.02,0.15,1.5]*1e-6;
sum(vThickness);
vNeffact=[3.1694,3.2525,3.31,3.364,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.364,3.31,3.2525,3.1694,3.1694];
vWidth=[totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,totalwid,3e-6];

%-----------% 八量子阱片%
vThickness=[1.5,0.025,0.025,0.025,0.025,tqw,tbarrier,tqw,tbarrier,tqw,tbarrier,tqw,tbarrier,tqw,tbarrier,tqw,tbarrier,tqw,tbarrier,tqw,0.025,0.025,0.025,0.15,1.5]*1e-6;
vNeffact=[3.1694,3.2525,3.2829,3.31,3.364,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.37,3.59,3.364,3.31,3.2525,3.1694,3.1694];
kk=length(vNeffact);
vWidth=[];
for i=1:(kk-2)
    vWidth=[vWidth,totalwid];
end
vWidth=[vWidth,3e-6,3e-6];

[eeps,nx,ny,mStEnd] = GenerateEps(vThickness,vNeffact,vWidth,ncld,totalwid,dx,dy);

lambda = 1.55e-6;
guess = 3.21;
% caculate the mode
boundary = 'oooo';  % 参svbuiltmx里面的说明
field = 'ey';
A = svbuildmtx (lambda,dx,dy,eeps,boundary,field);
bIsDrawing = 1;
nmodes=2;
[g,neff] = sveigenmodes(A,guess,nmodes,lambda,nx,ny,[],[],bIsDrawing);
neff



