function [alm] = scalarPlaneWaveCoef(L,thetak,phik)
% Scalar plane wave coefficients for fully normalized scalar
% wave functions
%
%   a_{lm} = 4\pi (i)^l Y_{l,m}^*(\theta_k,\phi_k)
%
% L:        maximum harmonic l = 0:L, m = -l:l
% thetak:   plane wave direction in theta (radians)
% phik:     plane wave direction in phi (radians)
%
% alm:      scalar plane wave coefficients linearly indexed
%
% Dependencies: lmtable, sphericalY

ylm = sphericalY(L,thetak,phik,'mono');
tab = lmtable(L,'mono');
l = tab(:,1);
alm = (4*pi)*(1i).^l.*conj(ylm(:));


