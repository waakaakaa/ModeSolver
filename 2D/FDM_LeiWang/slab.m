function [x,y,xc,yc,nx,ny,eps] = slab(n1,n2,n3,h1,h2,h3,h4,x1,x2,x3,wide,dx,dy);

% THis function creats an indes mesh for the finite-difference mode solver. 
% The function will accomodate a generalized three layer slab waveguide structure.

% Usage:
% [x,y,xc,yc,nx,ny,eps] = slab(n1,n2,n3,h1,h2,h3,wide,dx,dy);

% Input:
% n1 = index of refraction for lower layer
% n2 = index of refraction for the guide layer
% n3 = index of refraction for top layer
% h1 = height of lower layer
% h2 = height fo guide region
% h3 = height of top layer
% wide = the width of the slab in the calculation region 
% dx = horisontal grid spacing
% dy = vertical grid spacing

% Output:
% x,y = vectors specifying mesh coordinates
% xc,yc = vector specifying grid-center coordinates
% nx,ny = size of index mesh
% eps --- index mesh (n^2)


ih1 = round(h1/dy);%Àƒ…·ŒÂ»Î
ih2 = round(h2/dy);
ih3 = round(h3/dy);
ih4 = round(h4/dy);
ix1 = round(x1/dx);
ix2 = round(x2/dx);
ix3 = round(x3/dx);
iwide = round(wide/dx);
nx = ix1+ix2+ix3+1;
%nx = iwide+1;
ny = ih1+ih2+ih3+ih4+1;

xc = (1:(nx-1))'*dx-dx/2; %at the points 
yc = (1:(ny-1))'*dy-dy/2; 

x = (0:(nx-1))'*dx;%at the lapace
y = (0:(ny-1))'*dy;

eps = zeros(nx-1,ny-1);
% write N into the matrix%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iy=1;
ix=1;
for i = 1:ih1,
    eps(:,iy) = n1^2*ones(nx-1,1);
    iy = iy+1;
end

for i = 1:ih2,
    eps(:,iy) = n2^2*ones(nx-1,1);
    iy = iy+1;
end
uu=iy;
for i = 1:ih3,
    eps(:,iy) = n2^2*ones(nx-1,1);
    iy = iy+1;
end
% uu=iy-1;
for i = 1:ih4,
    eps(:,iy) = n3^2*ones(nx-1,1);
    iy = iy+1;
end

for i = 1:ih3,
    ix = 1;
    for j = 1:ix1,
        eps(ix,uu) = n3^2;
        ix = ix+1;
    end
    for j = 1:ix2,
        
        ix = ix+1;
    end
    for j = 1:ix3,
        eps(ix,uu) = n3^2;
        ix = ix+1;
    end
    uu = uu+1;
end
nx = length(xc);
ny = length(yc);