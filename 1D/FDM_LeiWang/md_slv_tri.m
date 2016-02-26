%3层波导的模式解
Pi = 3.1416;
epsl1 = 1;
epsl2 = -1;
L =0.05e-6;
lamda = 1.55e-6;
k0 = 2 * Pi/lamda;
I =sqrt(-1);
options=optimset('Display','off','TolX',1e-14,'LargeScale','off');
x0 = [sqrt(epsl1)*2;1];
[x,fval] = fsolve(@tri_lay_func,x0,options,epsl1,epsl2,L,lamda,k0)

k = k0*real(x(1))+I*real(x(2));
   k1 = sqrt(k0^2 * epsl1 - k.^2);
   if(real(k1)<0)
       k1 = -k1;
   end;
   k2 = sqrt(k0^2 * epsl2 - k.^2);
   if(imag(k2)>0)
       k2 = -k2;
   end;
b = 1;
a = (k1*epsl2 - k2*epsl1)/(k1*epsl2 + k2*epsl1);
c =a + b;
d = exp(-I*k1*L)*a+exp(I*k1*L)*b;

W = L*3;

dx = 0.001e-6;

n = round(W/dx);
Hy = zeros(1,n);
n1 = round((W-L)/2/dx);
n2 = round(((W-L)/2+L)/dx);

Hy(1:n1) = c*exp(I*k2*((1:n1)*dx - n1*dx));
Hy(n1+1:n2) = a*exp(-I*k1*((n1+1:n2)*dx - n1*dx))+b*exp(I*k1*((n1+1:n2)*dx - n1*dx));
Hy(n2+1:n) = d*exp(-I*k2*((n2+1:n)*dx - n1*dx - L));
plot(abs(Hy).^2);


