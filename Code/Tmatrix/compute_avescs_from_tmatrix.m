function [avescs] = compute_avescs_from_tmatrix(Tmm,Tmn,Tnm,Tnn,k)
% Compute polarization and orientation averaged scattering cross section from T-matrix
%
% Tmm,Tmn,Tnm,Tnn:    T-matrix block matrices
%                     Size [NxN], N = L^2 + 2L
% k:                  Background wavenumber
%
% avescs:             polarization and orientation averaged scattering cross section 

avescs = (2*pi/(abs(k)^2))*sum(abs(Tmm(:)).^2 + abs(Tmn(:)).^2 ...
    + abs(Tnm(:)).^2 + abs(Tnn(:)).^2);


