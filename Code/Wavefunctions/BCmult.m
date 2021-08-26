function [Fth, Fphi] = BCmult(Bth,Bphi,Cth,Cphi,blm,clm)
% Vector spherical function from vector spherical harmonics
%
% blm, clm:                 expansion coefficients
% Bth, Bphi, Cth, Cphi:     outputs from BC
%
% Fth, Fphi:                theta/phi field components

blm = blm(:);
clm = clm(:);
Fth = Bth*blm + Cth*clm;
Fphi = Bphi*blm + Cphi*clm;
