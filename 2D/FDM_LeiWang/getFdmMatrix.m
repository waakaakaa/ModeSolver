function [A] = getFdmMatrix (lambda,dx,dy,eeps,field)

[nx,ny] = size(eeps);

dyn = dy*ones(1,nx*ny);
dys = dy*ones(1,nx*ny);

eeps = [eeps(:,1),eeps,eeps(:,ny)];
eeps = [eeps(1,:);eeps;eeps(nx,:)];

k = 2*pi/lambda;

en = ones(1,nx*ny);
es = ones(1,nx*ny);
ee = ones(1,nx*ny);
ew = ones(1,nx*ny);
ep = ones(1,nx*ny);

en(:) = eeps(2:nx+1,3:ny+2);
es(:) = eeps(2:nx+1,1:ny);
ee(:) = eeps(3:nx+2,2:ny+1);
ew(:) = eeps(1:nx,  2:ny+1);
ep(:) = eeps(2:nx+1,2:ny+1);

switch lower(field)
    case 'te'
        alphaw = 4*(1+ep./ee)./(3+2*ep./ee+2*ep./ew+ep.^2./(ew.*ee));
        alphap = 2*(ep./ee+ep./ew+2*ep.^2./(ee.*ew))./(3+2*ep./ee+2*ep./ew+ep.^2./(ew.*ee));
        alphae = 4*(1+ep./ew)./(3+2*ep./ee+2*ep./ew+ep.^2./(ew.*ee));
        
        an = 2./(dyn.*(dyn+dys));
        as = 2./(dys.*(dyn+dys));
        ae = alphae/dx^2;
        aw = alphaw/dx^2;
        ap = ep*k^2-2*alphap/dx^2-2./(dyn.*dys);
        
    case 'tm'
        garmas = 4*(1+ep./en)./(3+2*ep./en+2*ep./es+ep.^2./(en.*es));
        garmap = 2*(ep./en+ep./es+2*ep.^2./(en.*es))./(3+2*ep./en+2*ep./es+ep.^2./(en.*es));
        garman = 4*(1+ep./es)./(3+2*ep./en+2*ep./es+ep.^2./(en.*es));
        
        an = garman/dy^2;
        as = garmas/dy^2;
        ae = ones(1,nx*ny)/dx^2;
        aw = ones(1,nx*ny)/dx^2;
        ap = ep*k^2-2*garmap/dy^2-2*ones(1,nx*ny)/dx^2;
end

ii = zeros(nx,ny);
ii(:) = (1:nx*ny);
iall = zeros(1,nx*ny);
is = zeros(1,nx*(ny-1));
in = zeros(1,nx*(ny-1));
ie = zeros(1,(nx-1)*ny);
iw = zeros(1,(nx-1)*ny);

iall(:) = ii;
is(:) = ii(1:nx,1:(ny-1));
in(:) = ii(1:nx,2:ny);
iw(:) = ii(1:(nx-1),1:ny);
ie(:) = ii(2:nx,1:ny);

A = sparse ([iall,iw,ie,is,in],[iall,ie,iw,in,is],[ap(iall),ae(iw),aw(ie),an(is),as(in)]);