clear
clc

% barrier
b1 = 0.01e-6;
n1 = 3.3647;

% well
b2 = 0.0055e-6;
n2 = 3.5241;

% number
n = 6;

nTE = sqrt( (n*b1*n1^2 + (n-1)*b2*n2^2) / (n*b1 + (n-1)*b2) );
nTM = 1 / sqrt( (n*b1*n1^(-2) + (n-1)*b2*n2^(-2)) / (n*b1 + (n-1)*b2) );