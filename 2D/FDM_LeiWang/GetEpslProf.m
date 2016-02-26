function [mEpsl,Nx,Ny,vNact,mGridInfo] = GetEpslProf(dx,dy1,dy2,dy3,Wdt,Wdtot,Dept)
%input

% dx --- mesh spacing in x axis
% dy1 --- mesh spacing in cladding and substrate sections
% dy2 ---- mesh spacing in GRINSCH sections
% dy3 ---- mesh spacing in QWs sections
% Wdt --- strip width
% Wdtot --- total width
% Dept --- etching depth of the grating section


%output
%mEpsl --- matrix of Epsl profile
%Nx,Ny --- grid numbers in x and y direction
%vNact --- vector denoting QW positions
% mGridInf --- the infomation of sections which have different grid size

Nx = round(Wdtot/dx);% the mesh number in x direction
Ny = 0;
mGridInfo = zeros(5,2);

%here we assume that the etch depth is deeper than the total thichness of p doped GRSINSCH layers
%--------------------the layered structure in x direction from the top to the bottom-------------------%
Thntot = 0;
I = sqrt(-1);
Thn_Cap = 0.2e-6; %P doped InGaAs cap layer
N_Cap = 3.5388 + 0.102*I;
ny = round(Thn_Cap/dy1);
mEpsl = N_Cap^2*ones(Nx,ny);
Ny = Ny + ny;
Thntot = Thntot + Thn_Cap;

Thn_Cld1 = 1.5e-6;%P doped InP cladding
N_Cld1 = 3.1694;
ny = round(Thn_Cld1/dy1);
mEpsl = [mEpsl,N_Cld1^2*ones(Nx,ny)];
Ny = Ny + ny;
N1 = Ny;
Thntot = Thntot + Thn_Cld1;

Thn_Es = 4e-9;%P doped InGaAs etch stop layer
N_Es = 3.389;

Thn_Cld2 = 0.15e-6;%P doped InP cladding
N_Cld2 = 3.1694;
ny = round(Thn_Cld2/dy1);
mEpsl = [mEpsl,N_Cld2^2*ones(Nx,ny)];
Ny = Ny + ny;
N2 = Ny+1;
Thntot = Thntot + Thn_Cld2;
Thn2 = Thntot;
mGridInfo(1,1) = Ny;
mGridInfo(1,2) = dy1;

Thn_GS1 = 0.025e-6;%P doped InGaAsP GRINSCH layer
N_GS1 = 3.21;
ny = round(Thn_GS1/dy2);
mEpsl = [mEpsl,N_GS1^2*ones(Nx,ny)];
Ny = Ny + ny;
Thntot = Thntot + Thn_GS1;

Thn_GS2 = 0.025e-6;%P doped InGaAsP GRINSCH layer
N_GS2 = 3.21;
ny = round(Thn_GS2/dy2);
mEpsl = [mEpsl,N_GS2^2*ones(Nx,ny)];
Ny = Ny + ny;
Thntot = Thntot + Thn_GS2;

Thn_GS3 = 0.025e-6;%P doped InGaAsP GRINSCH layer
N_GS3 = 3.21;
ny = round(Thn_GS3/dy2);
mEpsl = [mEpsl,N_GS3^2*ones(Nx,ny)];
Ny = Ny + ny;
Thntot = Thntot + Thn_GS3;
mGridInfo(2,1) = Ny;
mGridInfo(2,2) = dy2;

Thn_QW = 5.5e-9; %quantum well layers
N_QW = 3.59;
n_QW = 8; % the number of the QWs
vThnact = zeros(2*n_QW,1);

Thn_B = 10e-9
N_B = 3.37;% barrier layers
n_B = 7; % the number of the barriers
for i =1 : n_B
    vNact(i*2-1) = Ny+1;
    ny = round(Thn_QW/dy3);
    mEpsl = [mEpsl,N_QW^2*ones(Nx, ny)];
    Ny = Ny + ny;
    vNact(i*2) = Ny;
    
    ny = round(Thn_B/dy3);
    mEpsl = [mEpsl,N_B^2*ones(Nx,ny)];
    Ny = Ny + ny;
    Thntot = Thntot+Thn_QW+Thn_B;
end;
 i = i+1;
 vNact(i*2-1) = Ny+1;
 ny = round(Thn_QW/dy3);
 mEpsl = [mEpsl,N_QW^2*ones(Nx,ny)];
 Ny = Ny + ny;
 vNact(i*2) = Ny;
 mGridInfo(3,1) = Ny;
 mGridInfo(3,2) = dy3;
 

Thn_GS4 = 0.025e-6;%N doped InGaAsP GRINSCH layer
N_GS4 = 3.364;
ny = round(Thn_GS4/dy2);
mEpsl = [mEpsl,N_GS4^2*ones(Nx,ny)];
Ny = Ny + ny;

Thn_GS5 = 0.025e-6;%N doped InGaAsP GRINSCH layer
N_GS5 = 3.31;
ny = round(Thn_GS5/dy2);
mEpsl = [mEpsl,N_GS5^2*ones(Nx,ny)];
Ny = Ny + ny;

Thn_GS6 = 0.025e-6;%N doped InGaAsP GRINSCH layer
N_GS6 = 3.2829;
ny = round(Thn_GS6/dy2);
mEpsl = [mEpsl,N_GS6^2*ones(Nx,ny)];
Ny = Ny + ny;

Thn_GS7 = 0.025e-6;%N doped InGaAsP GRINSCH layer
N_GS7 = 3.2525;
ny = round(Thn_GS7/dy2);
mEpsl = [mEpsl,N_GS7^2*ones(Nx,ny)];
Ny = Ny + ny;
mGridInfo(4,1) = Ny;
mGridInfo(4,2) = dy2;

Thn_Sub = 2e-6;%N doped InP substrate layer(the refractive index of buffer is equal to that of the substrate, so it is negligible )
N_Sub = 3.1694;
ny = round(Thn_Sub/dy1);
mEpsl = [mEpsl,N_Sub^2*ones(Nx,ny)];
Ny = Ny + ny;
mGridInfo(5,1) = Ny;
mGridInfo(5,2) = dy1;

if(Dept > 0)
    N_Etch = round((Thn_GS1+Thn_GS2+Thn_GS3)/dy2+ (Dept-(Thn_GS1+Thn_GS2+Thn_GS3))/dy3);% mesh number of grating etching section
    mEpsl(:,N2:N2+N_Etch) = N_Cld1^2;
    for i =n_B+1:-1:1
       if(vNact(2*i)<N2+N_Etch)
               vNact = vNact(2*i+1:(n_B+1)*2);
           break;
       elseif (vNact(2*i-1)<N2+N_Etch)
           vNact(2*i-1)=N2+N_Etch;
           vNact = vNact(2*i-1:(n_B+1)*2);
           break;
       end;
    end;
end;

%-----------------------------------etching the bridge --------------------------%
Nbrg = round(Wdt/dx);
mEpsl(1:round((Nx-Nbrg)/2),1:N1) = 1;
mEpsl(round((Nx-Nbrg)/2)+Nbrg:Nx,1:N1) = 1;









