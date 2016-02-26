function [DEX,DEY,DHX,DHY] = yeeder(NGRID,RES,BC,kinc)
% YEEDER Construct Yee Grid Derivative Operators on a 2D Grid
%
% [DEX,DEY,DHX,DHY] = yeeder(NGRID,RES,BC,kinc);
%
% Note for normalized grid, use this function as follows:
%
% [DEX,DEY,DHX,DHY] = yeeder(NGRID,k0*RES,BC,kinc/k0);
%
% Input Arguments
% =================
% NGRID [Nx Ny] grid size
% RES [dx dy] grid resolution of the 1X grid
% BC [xbc ybc] boundary conditions
% -2: periodic (requires kinc)
% 0: Dirichlet
% kinc [kx ky] incident wave vector
% This argument is only needed for periodic boundaries.

Nx = NGRID(1);
Ny = NGRID(2);
dx = RES(1);
dy = RES(2);
xbc = BC(1);
ybc = BC(2);

ntotal = Nx*Ny;

if Nx == 1
    if xbc == -2
        B = 1i*kinc(1)*ones(1,ntotal)';
        DEX = 1/dx * sparse(spdiags(B, 0, ntotal,ntotal));
    else
        DEX = sparse(zeros(ntotal,ntotal));
    end
else
    B = 1/dx * [-ones(1,ntotal)' ones(1,ntotal)'];
    DEX = sparse(spdiags(B,[0 1],ntotal,ntotal));
    for i=Nx:Nx:Nx*Ny-1
        DEX(i,i+1) = 0;
    end
    if xbc == -2
        for i=Nx:Nx:ntotal
            DEX(i,i-Nx+1) = (1/dx)*exp(1i*kinc(1)*Nx*dx);
        end
    end
end
DHX = -DEX';

if Ny == 1
    if ybc == -2
        B = 1i*kinc(2)*ones(1,ntotal)';
        DEY = 1/dy * sparse(spdiags(B, 0, ntotal,ntotal));
    else
        DEY = sparse(zeros(ntotal,ntotal));
    end
else
    B = 1/dy * [-ones(1,ntotal)' ones(1,ntotal)'];
    DEY = sparse(spdiags(B,[0 Nx],ntotal,ntotal));
    if ybc == -2
        for i=1:Nx
            DEY(ntotal+i-Nx,i) = (1/dy)*exp(1i*kinc(2)*Ny*dy);
        end
    end
end
DHY = -DEY';
