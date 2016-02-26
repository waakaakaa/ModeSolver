%% ¿µÄþSMF-28  8.3/125¦Ìm,1.4681/1.4628
clear all
clc
close all

lambda0 = 0.55e-6;
k0 = 2*pi/lambda0;
ncore = 1.4681;
nclad = 1.4628;
nair = 1;
dcore = 8.3e-6;
dclad = 19e-6;

dx = 0.1e-6;
dy = 0.1e-6;
Nx = round(dclad/dx) +1;
Ny = Nx;

%%
Nx2 = 2*Nx;
dx2 = dx/2;
Ny2 = 2*Ny;
dy2 = dy/2;

center = round(dclad/2/dx2);
r1 = round(dcore/2/dx2);
r2 = center;

ER2 = nclad^2*ones(Nx2,Ny2);
for i = 1:Nx2
    for j = 1:Ny2
        if (i-center)^2 + (j-center)^2 <= r1^2
            ER2(i,j) = ncore^2;
        end
    end
end
UR2 = 1*ones(Nx2,Ny2);

subplot(2,2,1);
imagesc(UR2);colorbar;
subplot(2,2,2);
imagesc(ER2);colorbar;

URxx = UR2;
URyy = UR2;
URzz = UR2;
ERxx = ER2;
ERyy = ER2;
ERzz = ER2;

URxx = URxx(1:2:Nx2,2:2:Ny2);
URyy = URyy(2:2:Nx2,1:2:Ny2);
URzz = URzz(2:2:Nx2,2:2:Ny2);
ERxx = ERxx(2:2:Nx2,1:2:Ny2);
ERyy = ERyy(1:2:Nx2,2:2:Ny2);
ERzz = ERzz(1:2:Nx2,1:2:Ny2);

URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));

%%
[DEX,DEY,DHX,DHY] = yeeder([Nx Ny],[dx dy]*k0,[0 0]);

%%
P11 = DEX*ERzz^(-1)*DHY;
P12 = - (DEX*ERzz^(-1)*DHX + URyy);
P21 = DEY*ERzz^(-1)*DHY + URxx;
P22 = - DEY*ERzz^(-1)*DHX;
Q11 = DHX*URzz^(-1)*DEY;
Q12 = - (DHX*URzz^(-1)*DEX + ERyy);
Q21 = DHY*URzz^(-1)*DEY + ERxx;
Q22 = - DHY*URzz^(-1)*DEX;
OMEGA_square_PQ = [P11 P12;P21 P22]*[Q11 Q12;Q21 Q22];
NSOL = 1;
[V,D] = eigs(OMEGA_square_PQ, 2*NSOL, -ncore^2);
D = diag(D);
neff = - 1i*sqrt(D)
no = real(neff);
kappa = - imag(neff);
alpha = k0*kappa;
beta = k0*no;
VV1 = reshape(V(:,2*NSOL-1),Nx,Ny*2);
VV2 = reshape(V(:,2*NSOL),Nx,Ny*2);
subplot(2,2,3);
imagesc((0:dx:Nx*dx)*1e6,(0:dy:Ny*dy)*1e6,abs(VV2(:,1:Ny)));colorbar;
subplot(2,2,4);
mesh(abs(VV1(:,Ny+1:Ny*2)));