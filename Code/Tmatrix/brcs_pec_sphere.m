function [brcs] = brcs_pec_sphere(k,a)
% Backscatter Radar Cross Section of a PEC sphere
%
% k:    background wavenumber
% a:    sphere radius
%
% brcs:  Backscatter RCS of the PEC sphere
%
% Dependencies: sbesselj, sbesselh, sbesseljp2, sbesselhp2

ka = k*a;
l = 1;
err = 1;
old = 0;
while err > 1e-12;
    Tmm = sbesselj(l,ka)./sbesselh(l,ka);
    Tnn = sbesseljp2(l,ka)./sbesselhp2(l,ka);
    tmp = (-1).^l.*(2*l+1).*(Tmm - Tnn);
    brcs = old + tmp;
    err = abs(tmp)^2/abs(brcs)^2;
    old = brcs;
    l = l + 1;
end
brcs = pi/(abs(k)^2)*abs(brcs).^2;
