function [ v,d ] = getMode1D( eeps,k0,dx,guess,lambda0,mode )

n = length(eeps);
eps0 = 8.85e-12;
C0=2.997924e8;

epsl1 = eeps(1);
epsl2 = eeps(n);
nPMLwide = round(2e-6/dx);
amp = -11e3;

sigma_max1 = 0*amp/(k0*C0*eps0*epsl1);
sigma_max2 = 0*amp/(k0*C0*eps0*epsl2);

Sigma_XArray = zeros(1,n);
for i=1:nPMLwide
    Sigma_XArray(i) = (nPMLwide+0.5-(i-1))^2*sigma_max1/(nPMLwide)^2;
    Sigma_XArray(n-i) = (nPMLwide+0.5-(i-1))^2*sigma_max2/(nPMLwide)^2;
end;
if strcmp(mode, 'TE')
    n=n+1;
    Aip1 = ones(1,n-2)./((2.0 - 1i * Sigma_XArray(2:n-2+1) - 1i * Sigma_XArray(1:n-2)).*(1.0 - 1i * Sigma_XArray(2:n-2+1))/2*dx^2);
    Aim1 = ones(1,n-2)./((2.0 - 1i * Sigma_XArray(2:n-2+1) - 1i * Sigma_XArray(1:n-2)).*(1.0 - 1i * Sigma_XArray(1:n-2))/2*dx^2);
    Ai = k0^2*eeps(2:n-1)-2./((1.0 - 1i * Sigma_XArray(1:n-2)) .* ( 1.0 - 1i * Sigma_XArray(2:n-2+1))*dx^2);
    n=n-1;
elseif strcmp(mode, 'TM')
    Aip1 = 2*eeps(1:n-1)./(eeps(1:n-1) + eeps(2:n))./((2.0 - 1i * Sigma_XArray(2:n) - 1i * Sigma_XArray(1:n-1)).*(1.0 - 1i * Sigma_XArray(2:n))/2*dx^2);
    Aim1 = 2*eeps(2:n)./(eeps(1:n-1) + eeps(2:n))./((2.0 - 1i * Sigma_XArray(2:n) - 1i * Sigma_XArray(1:n-1)).*(1.0 - 1i * Sigma_XArray(1:n-1))/2*dx^2);
    Ai = 2*k0^2*eeps(1:n-1).*eeps(2:n)./(eeps(2:n)+eeps(1:n-1))-2./((1.0 - 1i * Sigma_XArray(1:n-1)) .* ( 1.0 - 1i * Sigma_XArray(2:n))*dx^2);
end;

n = n-1;
Aim1(1) = 0;
Aip1(n) = 0;

A = sparse([(1:n-1),(1:n),(2:n)], [(2:n),(1:n),(1:n-1)], [Aim1(2:n),Ai,Aip1(1:n-1)]);
shift = (2*pi*guess/lambda0)^2;
options.tol = 1e-24;
options.disp = 0;

[v,d] = eigs(A,1,shift,options);

end

