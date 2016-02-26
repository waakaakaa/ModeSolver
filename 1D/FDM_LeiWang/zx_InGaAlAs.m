clc
close all
clear all

lambda0 = 1.55e-6;
k0 = 2*pi/lambda0;
dx = 0.05e-9;
guess = 4; % the initial guess of the mode index

layers = [[1.2,3.17 ];...  %% substrate
    [0.01,  3.25    ];...  %% grin
    [0.06,  3.25    ];...  %% grin
    [0.06,  3.27    ];...  %% grin
    
    [0.01,   3.3647   ];...  %% barrier
    [0.0055, 3.5241   ];...  %% well
    [0.01,   3.3647   ];...  %% barrier
    [0.0055, 3.5241   ];...  %% well
    [0.01,   3.3647   ];...  %% barrier
    [0.0055, 3.5241   ];...  %% well
    [0.01,   3.3647   ];...  %% barrier
    [0.0055, 3.5241   ];...  %% well
    [0.01,   3.3647   ];...  %% barrier
    [0.0055, 3.5241   ];...  %% well
    [0.01,   3.3647   ];...  %% barrier
    [0.0055, 3.5241   ];...  %% well
    [0.01,   3.3647   ];...  %% barrier
    [0.0055, 3.5241   ];...  %% well
    [0.01,   3.3647   ];...  %% barrier
    [0.0055, 3.5241   ];...  %% well
    [0.01,   3.3647   ];...  %% barrier
    
    [0.06,   3.27   ];...    %% grin
    [0.06,   3.25   ];...    %% grin
    [0.06,   3.17   ];...    %% grin
    [0.01,   3.45   ];...    %% etch stop
    [1.5,    3.17   ]];...      %% cladding
    
[vThickness, vNeffAct] = getThickAndNFromLayers(layers);
[eeps] = getEpslFromThickAndN(vThickness,vNeffAct,dx);
[vEte, d] = getMode1D(eeps,k0,dx,guess,lambda0,'TE');
plotyy(sqrt(eeps),'r',abs(vEte.^2),'b');
neff = sqrt(d)/k0