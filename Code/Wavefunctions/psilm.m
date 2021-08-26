function [psi] = psilm(L,k,r,theta,phi,str)
% Scalar wave function:
%
%       psi_{lm}(kr) = h_l^(1)(kr)Y_l^m(theta,phi)
% or
%       Rg psi_{lm}(kr) = j_l(kr)Y_l^m(theta,phi)
%
% L:                Maximum harmonic degree L
% k:                Background wavenumber (real or complex)
% r,theta, phi:     Spherical coordinates, angles in radians
% str:              [optional] 'rg' for regular waves
%
% psilm:            Scalar wave function evaluated at (r,theta,phi)
%                   Dimensions: length(r(:)) x N
%                   N = L^2 + 2*L + 1
%
% Dependencies: lmtable, sphericalY, sbesselj, sbesselh

[l kr] = ndgrid(0:L,k*r(:));
if nargin == 6
    if ~strcmp(str,'rg')
        error('bad string')
    end
    bess = sbesselj(l,kr).';
else
    bess = sbesselh(l,kr).';
end
tab = lmtable(L,'mono');
psi = bess(:,tab(:,1)+1).*sphericalY(L,theta,phi,'mono'); 

