function [A,mdy] = svbuildmtx (lambda,dx,dy,eeps,boundary,field,mGridInfo);

% This function constructs the finite difference mattrix for the semivectorial mode solver.
% Using the index mesh specified in the matrix eeps(=n^2) and the boundary conditions specified in 'boundary',
% This function generates a sparse matrix A respresenting the differential operateor for the eigenmodes.

% Usage:
% A=svbuildmtx (lamda,dx,dy,eeps,boundary,field);

% Input:
% lambda -------- optical wavelength
% dx -------- horizental grid spacing
% dy -------- vertical grid spacing
% eeps ------- index mesh (=n^2(x,y))
% boundary ------ 4 letter string specifying boundary conditions to be
%    boundary(1) = North boundary condition
%    boundary(2) = South boundary condition
%    boundary(3) = East boundary condition
%    boundary(4) = West boundary condition
% The following boundary conditions are supported:
%  'A' ----- field is antisymmetric field outside of the boundary equal to the field at the boundary multiply -1
%  'S' ----- field is symmetric field outside of the boundary equal to the field at the boundary
%  'O' ------  field is zero immediately outside of the boundary
% Field ----- can be 'Ex', 'Ey', "Hx', 'Hy', or 'scale'
% mGridInf --- the infomation of sections which have different grid size

% Output:
% A ----- sparse matrix representing difference operator for the eigenvalue problem.
% mdy -----  the matrix denoting grid sizes in y direction


[nx,ny] = size(eeps);

% n = length(mGridInfo);
%
% mdy = mGridInfo(1,2)*ones(nx,mGridInfo(1,1));
% for i=2:n
%     mdy = [mdy,mGridInfo(i,2)*ones(nx,mGridInfo(i,1)-mGridInfo(i-1,1))];
% end;
%
dyn = dy*ones(1,nx*ny);
dys = dy*ones(1,nx*ny);
%
% dyn(:) = [mdy(:,2:ny),mdy(:,ny)];
% dys(:) = mdy;

% Now we pad eeps on all sides by one grid point
eeps = [eeps(:,1),eeps,eeps(:,ny)];
eeps = [eeps(1,:);eeps;eeps(nx,:)];

% Compute free-space wavevector
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
    case 'ex'% in y direction we introduce nonuniform grids
        alphaw = 4*(1+ep./ee)./(3+2*ep./ee+2*ep./ew+ep.^2./(ew.*ee));
        alphap = 2*(ep./ee+ep./ew+2*ep.^2./(ee.*ew))./(3+2*ep./ee+2*ep./ew+ep.^2./(ew.*ee));
        alphae = 4*(1+ep./ew)./(3+2*ep./ee+2*ep./ew+ep.^2./(ew.*ee));
        
        an = 2./(dyn.*(dyn+dys));
        as = 2./(dys.*(dyn+dys));
        ae = alphae/dx^2;
        aw = alphaw/dx^2;
        ap = ep*k^2-2*alphap/dx^2-2./(dyn.*dys);
        
    case 'ey'
        garmas = 4*(1+ep./en)./(3+2*ep./en+2*ep./es+ep.^2./(en.*es));
        garmap = 2*(ep./en+ep./es+2*ep.^2./(en.*es))./(3+2*ep./en+2*ep./es+ep.^2./(en.*es));
        garman = 4*(1+ep./es)./(3+2*ep./en+2*ep./es+ep.^2./(en.*es));
        
        an = garman/dy^2;
        as = garmas/dy^2;
        ae = ones(1,nx*ny)/dx^2;
        aw = ones(1,nx*ny)/dx^2;
        ap = ep*k^2-2*garmap/dy^2-2*ones(1,nx*ny)/dx^2;
        
    case 'hx'
        an = 2*ep./((ep+en)*dy^2);
        as = 2*ep./((ep+es)*dy^2);
        ae = ones(1,nx*ny)/dx^2;
        aw = ones(1,nx*ny)/dx^2;
        ap = ep*k^2-an-as-2*ones(1,nx*ny)/dx^2;
        
    case 'hy'
        an = ones(1,nx*ny)/dy^2;
        as = ones(1,nx*ny)/dy^2;
        ae = 2*ep./(((ep+ee)*dx^2)+eps);
        aw = 2*ep./(((ep+ew)*dx^2)+eps);
        ap = ep*k^2-ae-aw-2*ones(1,nx*ny)/dy^2;
        
    case 'scalar'
        an = ones(1,nx*ny)/dy^2;
        as = ones(1,nx*ny)/dy^2;
        ae = ones(1,nx*ny)/dx^2;
        aw = ones(1,nx*ny)/dx^2;
        ap = ep*k^2-2*ones(1,nx*ny)/dx^2-2*ones(1,nx*ny)/dy^2;
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


% Now we most account for the boundary conditions.

% North boundary
ib = zeros(1,nx);
b = boundary(1);
ib(:) = ii(1:nx,ny);
if (b == 's') sign = +1;
elseif (b == 'a') sign = -1;
elseif (b == 'o') sign = 0;
end
for i = ib,
    A(i,i) = A(i,i) + sign*an(i);
end

% South boundary
ib = zeros(1,nx);
b = boundary(2);
ib(:) = ii(1:nx,1);
if (b == 's') sign = +1;
elseif (b == 'a') sign = -1;
elseif (b == 'o') sign = 0;
end
for i = ib,
    A(i,i) = A(i,i) + sign*as(i);
end

% East boundary
ib = zeros(1,ny);
b = boundary(3);
ib(:) = ii(nx,1:ny);
if (b == 's') sign = +1;
elseif (b == 'a') sign = -1;
elseif (b == 'o') sign = 0;
end
for i = ib,
    A(i,i) = A(i,i) + sign*ae(i);
end

% West boundary
ib = zeros(1,ny);
b = boundary(4);
ib(:) = ii(1,1:ny);
if (b == 's') sign = +1;
elseif (b == 'a') sign = -1;
elseif (b == 'o') sign = 0;
end
for i = ib,
    A(i,i) = A(i,i) + sign*aw(i);
end