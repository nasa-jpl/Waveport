function [Tmm, Tnn] = tmatrixDielectricSphere(L,a,k1,k2,u1,u2)
% T-matrix of a dielectric sphere
%
% L:        Maximum harmonic degree L
% a:        radius of sphere (m)
% k1:       background wavenumber (1/m)
% k2:       sphere wavenumber (1/m)
% u1:       (optional) background relative permeability (must be included in k1)
% u2:       (optional) sphere relative permeability (must be included in k2)
%
% Tmm:      Nx1 diagonal elements of Tmatrix for Mlm, N = L^2 + 2*L	
% Tnn:      Nx1 diagonal elements of Tmatrix for Nlm
%
% Dependencies: lmtable, sbesselj, sbesselh, sbesslejp2, sbesselhp2

if nargin == 4
    u1 = 1;
    u2 = 1;
end
k1a = k1*a;
k2a = k2*a;
l = 1:L;
j1 = sbesselj(l,k1a);
j2 = sbesselj(l,k2a);
h1 = sbesselh(l,k1a);
j1p = sbesseljp2(l,k1a); 
j2p = sbesseljp2(l,k2a);
h1p = sbesselhp2(l,k1a);
N1 = u2*j2.*j1p - u1*j1.*j2p;
D1 = u1*h1.*j2p - u2*j2.*h1p;
N2 = u1*(k2a)^2*j2.*j1p - u2*(k1a)^2*j1.*j2p;
D2 = u2*(k1a)^2*h1.*j2p - u1*(k2a)^2*j2.*h1p;
Tmmtmp = N1./D1;
Tnntmp = N2./D2;
tab = lmtable(L);
Tmm = Tmmtmp(tab(:,1)).';
Tnn = Tnntmp(tab(:,1)).';