function [Er Eth Ephi] = MNmult(Mth,Mphi,Nr,Nth,Nphi,alm,blm)
% Spherical electric field from vector spherical wave functions
%
% alm, blm:                 expansion coefficients
% Mth,Mphi,Nr,Nth,Nphi:     outputs from MN
%
% Er,Eth,Ephi:              r/theta/phi field components

alm = alm(:);
blm = blm(:);
Er = Nr*blm;
Eth = Mth*alm + Nth*blm;
Ephi = Mphi*alm + Nphi*blm;


