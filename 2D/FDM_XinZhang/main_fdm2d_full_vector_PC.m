clear all
clc
close all

lambda0 = 1.55e-6;
k0 = 2*pi/lambda0;
nsi  = 1.5;
nair = 1.0;
dhole = 4.0e-6;
dh2h  = 4.5e-6;

dx = 0.05e-6;
dy = 0.05e-6;
Ny = round((2*dh2h+dhole)/dx);
Nx = round((2*dh2h*3^0.5/2+dhole)/dx);
%%
Nx2 = 2*Nx;
dx2 = dx/2;
Ny2 = 2*Ny;
dy2 = dy/2;

r = round(dhole/2/dx2);
c = [round(Nx2/2) round(dhole/2/dy2);...
    round(Nx2/2) round((dhole/2+2*dh2h)/dy2);...
    round(dhole/2/dx2) round(dhole/2/dx2)-round(dh2h/2/dx2);...
    round(dhole/2/dx2) round(dhole/2/dx2)+round(dh2h/2/dx2);...
    round(dhole/2/dx2) round(dhole/2/dx2)+round(dh2h/2/dx2)+round(dh2h/dx2);...
    round(dhole/2/dx2) round(dhole/2/dx2)+round(dh2h/2/dx2)+2*round(dh2h/dx2);...
    round(dhole/2/dx2)+round(dh2h*3^0.5/dx2) round(dhole/2/dx2)-round(dh2h/2/dx2);...
    round(dhole/2/dx2)+round(dh2h*3^0.5/dx2) round(dhole/2/dx2)+round(dh2h/2/dx2);...
    round(dhole/2/dx2)+round(dh2h*3^0.5/dx2) round(dhole/2/dx2)+round(dh2h/2/dx2)+round(dh2h/dx2);...
    round(dhole/2/dx2)+round(dh2h*3^0.5/dx2) round(dhole/2/dx2)+round(dh2h/2/dx2)+2*round(dh2h/dx2)];

ER2 = nsi^2*ones(Nx2,Ny2);
for k = 1:length(c)
    for i = 1:Nx2
        for j = 1:Ny2
            if (i-c(k,1))^2 + (j-c(k,2))^2 <= r^2
                ER2(i,j) = nair^2;
            end
        end
    end
end

UR2 = ones(Nx2,Ny2);

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
[V,D] = eigs(OMEGA_square_PQ, 2*NSOL, -(nsi-0.02)^2);
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