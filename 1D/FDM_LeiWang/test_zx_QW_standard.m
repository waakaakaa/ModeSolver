clc
close all
clear all

lambda0 = 1.55e-6;
k0 = 2*pi/lambda0;
dx = 0.5e-9;
boundary = 'oo';
nmodes = 1; % the number of modes to be caculated
guess = 4; % the initial guess of the mode index
mode = 'TE';

layers = [[1.5,    3.1694];...%% cladding
          [0.15,   3.1694];...%% etch stop
          [0.02,   3.2525];...%% grin
          [0.02,   3.31];...  %% grin
          [0.01,   3.364];... %% grin
          [0.01,   3.37];...  %% barrier
          [0.0055, 3.59];...  %% well
          [0.01,   3.37];...  %% barrier
          [0.0055, 3.59];...  %% well
          [0.01,   3.37];...  %% barrier
          [0.0055, 3.59];...  %% well
          [0.01,   3.37];...  %% barrier
          [0.0055, 3.59];...  %% well
          [0.01,   3.37];...  %% barrier
          [0.0055, 3.59];...  %% well
          [0.01,   3.37];...  %% barrier
          [0.01,   3.364];... %% grin
          [0.02,   3.31];...  %% grin
          [0.02,   3.2525];...%% grin
          [30,     3.1694]];  %% substrate

[vThickness, vNeffAct] = getThickAndNFromLayers(layers);
[eeps] = getEpslFromThickAndN(vThickness,vNeffAct,dx);
[vEte, d] = getMode1D(eeps,k0,dx,guess,lambda0,'TE');
