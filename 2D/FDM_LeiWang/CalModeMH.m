
% the perimittity profile, the mesh number in x and y direction of the modulation section
% vNact_y ---- the vector denoting the active layers(quantum wells)
% the perimittity profile, the mesh number in x and y direction of the modulation section
% vNact ---- the vector denoting the active layers(quantum wells)
[mEpsl_MH,Nx_M,Ny_M,vNact,mGridInfo] = GetEpslProf(dx,dy1,dy2,dy3,Wdt_M,Wdtot,0);


field = 'ex';
vboundary = 'oooo';
nmodes= 1;


vlambda = lambda_B+[-10:10]*1e-9;
nlam = length(vlambda);
vGama_MH = zeros(1,nlam);% the confinement factor of the modulation section 
vNeff0_MH = zeros(1,nlam);% the effective index without absorption of modulation section


for(i=-0:0)
    i
lambda = vlambda(i+11);
[A,mdy] = svbuildmtx (lambda,dx,dy,mEpsl_MH,vboundary,field,mGridInfo);% the difference 
guess = 3.3;
[mE_Prof_MH,vGama_MH(i+11),vNeff0_MH(i+11)] = sveigenmodes(A,guess,nmodes,lambda,Nx_M,Ny_M,vNact,mdy);
sFilename = GetFileName('E:/semiconductor laser/program/modulation/data/','mE_Prof_MH',Wdt_M,0,i+11);

save(sFilename,'mE_Prof_MH');
end;
sFilename = GetFileName('E:/semiconductor laser/program/modulation/data/','mdy',Wdt_M,0,0);
save(sFilename,'mdy');
save(GetFileName('E:/semiconductor laser/program/modulation/data/','vGama_MH',Wdt_M,0,0),'vGama_MH');
save(GetFileName('E:/semiconductor laser/program/modulation/data/','vNeff0_MH',Wdt_M,0,0),'vNeff0_MH');

