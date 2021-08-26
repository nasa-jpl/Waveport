function [alm, blm] = vectorPlaneWaveCoefZ(L,Ex,Ey)
% Multipole coefficients for z-propagating plane wave at the orgin
% for fully-normalized vector wave functions.
%
% L:            Maximum harmonic L
% Ex,Ey:        Electric field components
%
% alm, blm:     z-hat plane wave multipole coefficients, 
%               harmonics 1:L, all m, linearly indexed

tot = L^2 + 2*L;
alm = zeros(tot,1);
blm = zeros(tot,1);
for l=1:L,
    tmp = (1i)^(l+1)*sqrt(pi*(2*l+1));
    ind = l^2 + l + 1;
    alm(ind) = tmp*(Ex - 1i*Ey);
    blm(ind) = alm(ind);
    ind = l^2 + l - 1;
    alm(ind) = tmp*(Ex + 1i*Ey);
    blm(ind) = -alm(ind);
end

