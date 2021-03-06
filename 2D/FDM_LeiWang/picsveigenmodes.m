function [phi,neff]=sveigenmodes(A,guess,nmodes,lambda,nx,ny);
%function [neff]=sveigenmodes(A,guess,nmodes,lambda,nx,ny);

% This function calcualtes and formate the eigenmodes for the semivectorial finite difference mode solver.
% It uses the MATLAB function eigs to calculates a few eigenmodes of the matrix A.

% Usage:
% [phi,neff]=sveigenmodes(A,guess,nmodes,lambda,nx,ny);

% Input:
% A ----- sparse matrix containing the finite difference representation og the difference operator for the mode solver, generated by avbuildmtx.m
% guess --- scalar shift to apply when calculating the eigenvalues. This routine will return the eigenpairs which are closeat th this guess in magnitude.
% nmodes ---- number of modes to be calculated.
% lambda --- wavelength
% nx,ny ---- dimension of finite difference mesh

% Output:
% phi ---- three-dimensional vector containing the field for each calculated mode, phi(:,:,k)= kth eigenmode
% neff --- vector of modal effective indices, neff(k) = effective index of kth eigenmode


shift = (2*pi*guess/lambda)^2;
options.tol = 1e-8;
options.disp = 0;
[v,d] = eigs(A,speye(size(A)),nmodes,'lr',options);
%[v,d] = eigs(A);
neff = lambda*sqrt(diag(d))/(2*pi);
phi = v;
%temp=linspace(0,0,tte);
%for k = 1:nmodes
 %  phi(:,:,k) = v(:,:,k)/max(abs(v(:,:,k)));
  %  phi(:,:,k) = temp;
  % end;
for k=1:nmodes
    Hx = reshape(abs(v(1:nx*ny,k)),[nx ny]);
    figure;
    mesh((Hx));
    colorbar;
    
    Hy = reshape(abs(v(1+nx*ny:2*nx*ny,k)),[nx ny]);
    figure;
    mesh((Hy));
    colorbar;
    
    H=sqrt(Hx.^2+Hy.^2);
    figure;
    mesh((H));
    colorbar;
end;