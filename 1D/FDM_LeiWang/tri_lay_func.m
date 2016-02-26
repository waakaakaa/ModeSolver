function  y = tri_lay_func(x,epsl1,epsl2,L,lamda,k0);
   I =sqrt(-1);
   k = k0*real(x(1))+I*real(x(2));
   k1 = sqrt(k0^2 * epsl1 - k.^2);
   if(real(k1)<0)
       k1 = -k1;
   end;
   k2 = sqrt(k0^2 * epsl2 - k.^2);
   if(imag(k2)>0)
       k2 = -k2;
   end;
   y = ((k1*epsl2+k2*epsl1)/(k1*epsl2-k2*epsl1)).^2-exp(-2*I*k1*L);