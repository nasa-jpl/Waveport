function [Mth Mphi Nr Nth Nphi] = MN(L,k,r,theta,phi,rgstr,hatstr,normstr)
% Vector wave functions:
%
%       M_lm(kr) = h_l^(1)(kr)C_lm(theta,phi)
% and
%       N_lm(kr) = l(l+1)h_l^(1)(kr)/kr P_lm(theta,phi)
%                    + [kr h_l^(1)(kr)]'/kr B_lm(theta,phi)
% or
%
%       RgM_lm(kr) = j_l(kr)C_lm(theta,phi)
% and
%       RgN_lm(kr) = l(l+1)j_l(kr)/kr P_lm(theta,phi)
%                      + [kr j_l(kr)]'/kr B_lm(theta,phi)
%
% L:                Maximum harmonic degree L
% k:                Background wavenumber (real or complex)
% r,theta,phi:      Spherical coordinates, angles in radians
% str:              'rg' for regular waves
% hat:              'hat' for angular conjugation: \hatM and \hatN
% normstr:          'norm' for full, 1/sqrt(1(l+1)), normalization
%
% Mth,Mphi,         Vector wave function evaluated at (r,theta,phi)
% Nr,Nth,Nphi:      Dimensions: lenght(r(:)) x N
%                   N = L^2 + 2*L
%
% Dependencies: ind2lm, sphericalY, sbesselj, sbesselj2, 
%               sbesseljp, sbesselh, sbesselhp, BC

N = L^2 + 2*L;
[l m] = ind2lm((1:N)');
kr = k*r(:);
[l2 KR] = ndgrid(1:L,kr);
if nargin < 5
    error('not enough inputs')
end
if nargin >= 6 && strcmp(rgstr,'rg')
    bess1 = sbesselj(l2,KR).';      % j_l(kr), handles lim r->0 
    bess2 = sbesselj2(l2,KR).';     % j_l(kr)/kr, "
    bess3 = sbesseljp(l2,KR).';     % j'_l(kr), " 
else
    if ~isempty(rgstr)
        error('bad string')
    end
    bess1 = sbesselh(l2,KR).';          % h_l^(1)(kr)
    bess2 = (sbesselh(l2,KR)./KR).';    % h_l^(1)(kr)/kr
    bess3 = sbesselhp(l2,KR).';         % h'_l^(1)(kr)
end
bess4 = bess2 + bess3;
if nargin == 8 && strcmp(normstr,'norm')
    ylm = repmat(sqrt(l.*(l+1))',length(kr),1).*sphericalY(L,theta,phi);
    [Bth, Bphi, Cth, Cphi] = BC(L,theta,phi,'norm');
elseif nargin == 8 && ~isempty(normstr)
    error('bad string')
else
    ylm = repmat((l.*(l+1))',length(kr),1).*sphericalY(L,theta,phi);
    [Bth, Bphi, Cth, Cphi] = BC(L,theta,phi);
end
if nargin >= 7 && strcmp(hatstr,'hat')
    ylm = conj(ylm);
    Bth = conj(Bth);
    Bphi = conj(Bphi);
    Cth = conj(Cth);
    Cphi = conj(Cphi);
else
    if ~isempty(hatstr), 
        error('bad string')
    end
end
Mth = bess1(:,l).*Cth;        
Mphi = bess1(:,l).*Cphi;
Nr = bess2(:,l).*ylm;
Nth = bess4(:,l).*Bth;
Nphi  = bess4(:,l).*Bphi;