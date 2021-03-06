function [v,d] = Modesolver1D(epsl,k0,dx,nmodes,guess,lamda,boundary,Mode)%1-D mode solver for Hy
% [vAct dAct] = Modesolver1D(eepsAct,k0,dx,nmodes,guess,lambda,boundary,'TE');
    PI = 3.1415926;
    n = length(epsl);
    I =sqrt(-1);
    
    eps0 = 8.85e-12;
    C0=2.997924e+8;
    %in x direction
    epsl1 = epsl(1);
    epsl2 = epsl(n);
    nPMLwide = round(2e-6/dx);
    reflectivity = 1e-5;
    amp = -11e3;

    sigma_max1 = 0*amp/(k0*C0*eps0*epsl1);
    sigma_max2 = 0*amp/(k0*C0*eps0*epsl2);
  
  %������sigma_max�Ǹ��ģ�����%
  
    Sigma_XArray = zeros(1,n);
	for (i=1:nPMLwide)
		Sigma_XArray(i) = (nPMLwide+0.5-(i-1))^2*sigma_max1/(nPMLwide)^2;
		Sigma_XArray(n-i) = (nPMLwide+0.5-(i-1))^2*sigma_max2/(nPMLwide)^2;
	end;
    if(Mode == 'TE')
        n=n+1;
        Aip1 = ones(1,n-2)./((2.0 - I * Sigma_XArray(2:n-2+1) - I * Sigma_XArray(1:n-2)).*(1.0 - I * Sigma_XArray(2:n-2+1))/2*dx^2);%;epsl(1:n-2)./(epsl(1:n-2) + epsl(3:n))./((2.0 - I * Sigma_XArray(2:n-2+1) - I * Sigma_XArray(1:n-2)).*(1.0 - I * Sigma_XArray(2:n-2+1))/2*dx^2);%-c1
        Aim1 = ones(1,n-2)./((2.0 - I * Sigma_XArray(2:n-2+1) - I * Sigma_XArray(1:n-2)).*(1.0 - I * Sigma_XArray(1:n-2))/2*dx^2);%epsl(3:n)./(epsl(1:n-2) + epsl(3:n))./((2.0 - I * Sigma_XArray(2:n-2+1) - I * Sigma_XArray(1:n-2)).*(1.0 - I * Sigma_XArray(1:n-2))/2*dx^2);%c2;
        Ai = k0^2*epsl(2:n-1)-2./((1.0 - I * Sigma_XArray(1:n-2)) .* ( 1.0 - I * Sigma_XArray(2:n-2+1))*dx^2);%k0^2*(epsl(1:n-2).*epsl(3:n))./(epsl(1:n-2) + epsl(3:n))-1./((1.0 - I * Sigma_XArray(1:n-2)) .* ( 1.0 - I * Sigma_XArray(2:n-2+1))*dx^2);%c1.*d2 - c2.*d1;
        n=n-1;
    elseif(Mode == 'TM')
       Aip1 = 2*epsl(1:n-1)./(epsl(1:n-1) + epsl(2:n))./((2.0 - I * Sigma_XArray(2:n) - I * Sigma_XArray(1:n-1)).*(1.0 - I * Sigma_XArray(2:n))/2*dx^2);
       Aim1 = 2*epsl(2:n)./(epsl(1:n-1) + epsl(2:n))./((2.0 - I * Sigma_XArray(2:n) - I * Sigma_XArray(1:n-1)).*(1.0 - I * Sigma_XArray(1:n-1))/2*dx^2);     
       Ai = 2*k0^2*epsl(1:n-1).*epsl(2:n)./(epsl(2:n)+epsl(1:n-1))-2./((1.0 - I * Sigma_XArray(1:n-1)) .* ( 1.0 - I * Sigma_XArray(2:n))*dx^2);%k0^2*(epsl(1:n-2).*epsl(3:n))./(epsl(1:n-2) + epsl(3:n))-1./((1.0 - I * Sigma_XArray(1:n-2)) .* ( 1.0 - I * Sigma_XArray(2:n-2+1))*dx^2);%c1.*d2 - c2.*d1;
    end;
%      lm = 1/2;
%     for(i=1:n-2)
%         if(epsl(i+2) == epsl(i))
%            Aip1(i) = 1/dx^2;
%            Aim1(i) = 1/dx^2;
%            Ai(i) = -2/dx^2 + k0^2*epsl(i+1);
%        else
%          if(epsl(i+1)==epsl(i+2))
%             A =[-(lm+1)^2*epsl(i), -(lm-1)^2*epsl(i+2), lm^2*epsl(i),      epsl(i);
%             -(lm+1)^2,           -(lm-1)^2,         lm^2,              0;
%             -(lm)^2*epsl(i), -(lm-2)^2*epsl(i+2),  (lm-1)^2*epsl(i+2), epsl(i+2);
%             -(lm)^2,           -(lm-2)^2,         (lm-1)^2,            0;
%            ];
%          else
%             A =[-(lm-1)^2*epsl(i+2), -(lm+1)^2*epsl(i), lm^2*epsl(i),      epsl(i);
%             -(lm-1)^2,           -(lm+1)^2,         lm^2,              0;
%             -(lm-2)^2*epsl(i+2), -(lm)^2*epsl(i),  (lm-1)^2*epsl(i+2), epsl(i+2);
%              -(lm-2)^2,           -(lm)^2,         (lm-1)^2,            0;
%            ]; 
%           end;
%        
%         B = [0;-1;0;-1];
%         C = inv(A)*B;
%         Aip1(i) = 2*C(2)/dx^2;
%         Aim1(i) = 2*C(1)/dx^2;
%         Ai(i) = -2*C(3)/dx^2 + k0^2*epsl(i+1)*C(4);

    n = n-1;
    Aim1(1) = 0;
    Aip1(n) = 0;
    if(boundary(1) == 'o')%zero boundary condition
        
    else
        if(boundary(1) == 's')%symmetry boundary condition
            Ai(1) = Ai(1)+epsl(2)./(epsl(1) + epsl(2))/dx^2;
        else    %asymmetry  boundary condition
             Ai(1) = Ai(1)-epsl(2)./(epsl(1) + epsl(2))/dx^2;
        end;
    end;
    
    if(boundary(2) == 'o')%zero boundary condition
        
    else
        if(boundary(2) == 's')%symmetry boundary condition
            Ai(n) = Ai(n)+epsl(n)./(epsl(n) + epsl(n+1))/dx^2;
        else    %asymmetry  boundary condition
             Ai(n) = Ai(n);
        end;
    end;
    
    iall = [1:n];
    A = sparse([iall(1:n-1),iall,iall(2:n)],[iall(2:n),iall,iall(1:n-1)],[Aim1(2:n),Ai,Aip1(1:n-1)]);
    shift = (2*PI*guess/lamda)^2;
    options.tol = 1e-24;
    options.disp = 0;
    
    [v,d,flag] = eigs(A,nmodes,shift,options);