clc
clear

dx       = 0.02e-6;
dy       = 0.002e-6;
totalwid = 12;
wgwid    = 3;
lambda   = 1.55e-6;

layers = [[2,3.1694, totalwid];...  %% substrate
    [0.025,  3.2525, totalwid];...  %% grin
    [0.025,  3.2829, totalwid];...  %% grin
    [0.025,  3.31,   totalwid];...  %% grin
    [0.025,  3.364,  totalwid];...  %% grin
    
    [0.0055, 3.59,   totalwid];...  %% well
    [0.01,   3.37,   totalwid];...  %% barrier
    [0.0055, 3.59,   totalwid];...  %% well
    [0.01,   3.37,   totalwid];...  %% barrier
    [0.0055, 3.59,   totalwid];...  %% well
    [0.01,   3.37,   totalwid];...  %% barrier
    [0.0055, 3.59,   totalwid];...  %% well
    [0.01,   3.37,   totalwid];...  %% barrier
    [0.0055, 3.59,   totalwid];...  %% well
    [0.01,   3.37,   totalwid];...  %% barrier
    [0.0055, 3.59,   totalwid];...  %% well
    [0.01,   3.37,   totalwid];...  %% barrier
    [0.0055, 3.59,   totalwid];...  %% well
    [0.01,   3.37,   totalwid];...  %% barrier
    [0.0055, 3.59,   totalwid];...  %% well
    
    [0.025,  3.364,  totalwid];...  %% grin
    [0.025,  3.31,   totalwid];...  %% grin
    [0.025,  3.2525, totalwid];...  %% grin
    [0.15,   3.1694, wgwid];...     %% grin
    [1.5,    3.1694, wgwid];...     %% cladding
    [0.5,    1.0,    totalwid]];... %% air
    
[vThickness, vNeffact, vWidth] = getThickAndNAndWidFromLayers(layers);

[eeps,nx,ny] = getEpsl(vThickness,vNeffact,vWidth,totalwid*1e-6,dx,dy);

A = getFdmMatrix (lambda,dx,dy,eeps,'TE');

[mE,neff_te] = getModes2D(A,max(vNeffact),lambda,nx,ny);

imagesc(linspace(-totalwid,totalwid,nx),[0,sum(vThickness)*1e6],mE);set(gca,'ydir','normal');colorbar;

% optical confinement factor
layerNumbers = round(vThickness/dy);
Gama = zeros(1,length(layerNumbers));
for k = 1:length(layerNumbers)-1
    istart = k;
    iend   = k+1;
    nstart = sum(layerNumbers(1:(istart-1))) + 1;
    nend   = sum(layerNumbers(1:(iend  -1)));
    Gama(1,k) = sum(sum(abs(mE(nstart:nend,:).^2)))/sum(sum(abs(mE(:,:).^2)));
end
sum(Gama(1,6:20))