function [s,neff] = getModes2D(A,guess,lambda,nx,ny)

shift = (2*pi*guess/lambda)^2;
options.tol = 1e-8;
options.disp = 0;
[v,d] = eigs(A,speye(size(A)),1,shift,options);
neff = lambda*sqrt(d)/(2*pi);
s = reshape(v(:,1),[nx ny]);
s = s/max(max(s));
s = s';