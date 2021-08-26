function [Tmm, Tnn] = tmatrixPECSphere(L,a,k)
% T-matrix of a PEC sphere
%
% L:        Maximum harmonic degree L
% a:        radius of sphere (m)
% k:        background wavenumber (1/m)
% Tmm:      Nx1 diagonal elements of Tmatrix for Mlm, N = L^2 + 2*L	
% Tnn:      Nx1 diagonal elements of Tmatrix for Nlm
%
% Dependencies: lmtable, sbesselj, sbesselh, sbesslejp2, sbesselhp2

ka = k*a;
l=1:L;
Tmmtmp = -sbesselj(l,ka)./sbesselh(l,ka);
Tnntmp = -sbesseljp2(l,ka)./sbesselhp2(l,ka);
tab = lmtable(L);
Tmm = Tmmtmp(tab(:,1)).';
Tnn = Tnntmp(tab(:,1)).';


