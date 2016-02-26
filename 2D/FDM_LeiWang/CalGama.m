function Gama = CalGama(mIn_Prof,dx,dy,Nx,Ny,vNact,mdy)
% the function returns the confinement factor of the mode intensity profile
% described by matrix mIn_Prof

Ptot = 0;
Pact = 0;
n = length(vNact);

dyn = [mdy(:,2:Ny),mdy(:,Ny)];
dys = mdy;
temp = mIn_Prof.*(dyn+dys)/2;
Ptot = sum(sum(temp));
for j = 1:2:n-1
    Pact = Pact + sum(sum(temp(:,vNact(j):vNact(j+1))));
end;
% for i=1:Nx
%    Ini = interp1([0:Ny-1],mIn_Prof(i,:),[0:1/ratio:Ny-1]); %insert values in y axis
%
%    Ptot = Ptot + sum(Ini);
% end;
Gama = Pact/Ptot;
