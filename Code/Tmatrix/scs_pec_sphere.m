function [scs] = scs_pec_sphere(k,a)
% Scattering Cross Section of a PEC Sphere
%
% k:	background wavenumber
% a:	sphere radius
%
% scs:  scattering cross section of the PEC sphere
%
% Dependencies: sbesselj, sbesselh, sbesseljp2, sbesselhp2

x = k*a;
l = 1;
err = 1;
old = 0;
scs = 0;
while err > eps;
    j1 = sbesselj(l,x);
    j1p = sbesseljp2(l,x);
    h1 = sbesselh(l,x);
    h1p = sbesselhp2(l,x);
    t1 = abs(j1./h1).^2;
    t2 = abs(j1p./h1p).^2;
    tmp = (2*l+1).*(t1 + t2);
    scs = old + tmp;
    err = abs(tmp)^2/abs(scs)^2;
    old = scs;
    l = l + 1;
end
scs = 2*pi/(abs(k)^2)*scs;