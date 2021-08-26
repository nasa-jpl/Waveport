function [brcs] = brcs_dielectric_sphere(k1,k2,a)
% Backscatter Radar Cross Section of a Dielectric Sphere
%
% k1:   background wavenumber
% k2:   sphere wavenumber   
% a:    sphere radius
%
% brcs:  Backscatter RCS of the dielectric sphere
%
% Dependencies: sbesselj, sbesselh, sbesseljp2, sbesselhp2

x1 = k1*a;
x2 = k2*a;
l = 1;
err = 1;
old = 0;
brcs = 0;
while err > eps;
    j2 = sbesselj(l,x2);
    j2p = sbesseljp2(l,x2);
    h1 = sbesselh(l,x1);
    h1p = sbesselhp2(l,x1);
    num = j2.*j2p.*(x2.^2 - x1.^2)./x1;
    t1 = h1.*j2p;
    t2 = j2.*h1p;
    den = (t1 - t2).*(x1.^2.*t1 - x2.^2.*t2);
    tmp = (-1).^l.*(2*l+1).*(num./den);
    brcs = old + tmp;
    err = abs(tmp)^2/abs(brcs)^2;
    old = brcs;
    l = l + 1;
end
brcs = pi/(abs(k1)^2)*abs(brcs).^2;
