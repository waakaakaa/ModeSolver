clc
clear

dx     = 0.02e-6;
dy     = 0.02e-6;
lambda = 1.15e-6;

figure(8);hold on;box on;
xlabel('t (um)');ylabel('normalized propagation constants');
xlim([0 1]);ylim([0.25 0.4]);

for t = [0.1 0.3 0.5 0.7 0.9]
    layers = getLayers(t);
    [vThickness, vNeffact, vWidth] = getThickAndNAndWidFromLayers(layers);
    [eeps,nx,ny] = getEpsl(vThickness,vNeffact,vWidth,8e-6,dx,dy);
    
    A = getFdmMatrix (lambda,dx,dy,eeps,'TE');
    [s,neff] = getModes2D(A,max(vNeffact),lambda,nx,ny);
    b = (neff^2 - 3.4^2)/(3.44^2-3.4^2);
    plot(t,b,'ro');
    text(t,b,['   ',num2str(b)])
    
    A = getFdmMatrix (lambda,dx,dy,eeps,'TM');
    [s,neff] = getModes2D(A,max(vNeffact),lambda,nx,ny);
    b = (neff^2 - 3.4^2)/(3.44^2-3.4^2);
    plot(t,b,'bo');
    text(t,b,['   ',num2str(b)])
end

hold off;