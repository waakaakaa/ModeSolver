clc
clear

dx       = 0.02e-6;
dy       = 0.002e-6;
totalwid = 12;
wgwid    = 3;
lambda   = 1.55e-6;

layers = [[2.3,3.17, totalwid];...  %% substrate
    [0.01,  3.25,    totalwid];...  %% grin
    [0.06,  3.25,    totalwid];...  %% grin
    [0.06,  3.27,    totalwid];...  %% grin
    
    [0.01,   3.3647,   totalwid];...  %% barrier
    [0.0055, 3.5241,   totalwid];...  %% well
    [0.01,   3.3647,   totalwid];...  %% barrier
    [0.0055, 3.5241,   totalwid];...  %% well
    [0.01,   3.3647,   totalwid];...  %% barrier
    [0.0055, 3.5241,   totalwid];...  %% well
    [0.01,   3.3647,   totalwid];...  %% barrier
    [0.0055, 3.5241,   totalwid];...  %% well
    [0.01,   3.3647,   totalwid];...  %% barrier
    [0.0055, 3.5241,   totalwid];...  %% well
    [0.01,   3.3647,   totalwid];...  %% barrier
    [0.0055, 3.5241,   totalwid];...  %% well
    [0.01,   3.3647,   totalwid];...  %% barrier
    [0.0055, 3.5241,   totalwid];...  %% well
    [0.01,   3.3647,   totalwid];...  %% barrier
    [0.0055, 3.5241,   totalwid];...  %% well
    [0.01,   3.3647,   totalwid];...  %% barrier
    
    [0.06,   3.27,   totalwid];...    %% grin
    [0.06,   3.25,   totalwid];...    %% grin
    [0.06,   3.17,   totalwid];...    %% grin
    [0.01,   3.45,   totalwid];...    %% etch stop
    [1.5,    3.17,   wgwid]];...      %% cladding
    
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
sum(Gama(1,5:21))