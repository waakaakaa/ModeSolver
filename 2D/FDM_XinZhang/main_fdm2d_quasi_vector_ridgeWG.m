clear all
clc
close all

slab_height = 0.5e-6;

lambda0 = 1.15e-6;
k0 = 2*pi/lambda0;
ncore = 3.44;
nsub = 3.40;
nair = 1.0;
ridge_width = 3e-6;
ridge_height = 1e-6 - slab_height;
slab_width = 8e-6;
sub_height = 3e-6;
air_height = 0.5e-6;

dx = 0.02e-6;
dy = 0.02e-6;
Nx = round((air_height + ridge_height + slab_height + sub_height)/dx) +1;
Ny = round(slab_width/dy) +1;
%%
Nx2 = 2*Nx;
dx2 = dx/2;
Ny2 = 2*Ny;
dy2 = dy/2;

nx1 = 1 + round(air_height/dx2);
nx2 = nx1 + round(ridge_height/dx2) - 1;
nx3 = nx2 + 1;
nx4 = nx3 + round(slab_height/dx2) - 1;
nx5 = nx4 + 1;

ny1 = 1 + round((slab_width-ridge_width)/2/dy2);
ny2 = ny1 + round(ridge_width/dy2) - 1;

ER2 = nsub^2*ones(Nx2,Ny2);
ER2(1:nx2,:) = nair^2;
ER2(nx3:nx4,:) = ncore^2;
ER2(nx1:nx2,ny1:ny2) = ncore^2;

UR2 = 1*ones(Nx2,Ny2);

subplot(4,1,1);
imagesc(UR2);colorbar;
subplot(4,1,2);
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
OMEGA_xx_square = DEX*ERzz^(-1)*DHY*DHX*URzz^(-1)*DEY ...
    - (DEX*ERzz^(-1)*DHX + URyy) ...
    * (DHY*URzz^(-1)*DEY + ERxx);
NSOL = 1;
[V,D] = eigs(OMEGA_xx_square, NSOL, -ncore^2);
D = diag(D);
neff = - 1i*sqrt(D)
no = real(neff);
kappa = - imag(neff);
alpha = k0*kappa;
beta = k0*no;
VV = reshape(V(:,NSOL),Nx,Ny);
subplot(4,1,3);
imagesc((0:dx:Nx*dx)*1e6,(0:dy:Ny*dy)*1e6,abs(VV));colorbar;
%%
OMEGA_yy_square = DEY*ERzz^(-1)*DHX*DHY*URzz^(-1)*DEX ...
    - (DEY*ERzz^(-1)*DHY + URxx) ...
    * (DHX*URzz^(-1)*DEX + ERyy);
NSOL = 1;
[V,D] = eigs(OMEGA_yy_square, NSOL, -ncore^2);
D = diag(D);
neff = - 1i*sqrt(D)
no = real(neff);
kappa = - imag(neff);
alpha = k0*kappa;
beta = k0*no;
VV = reshape(V(:,NSOL),Nx,Ny);
subplot(4,1,4);
imagesc((0:dx:Nx*dx)*1e6,(0:dy:Ny*dy)*1e6,abs(VV));colorbar;