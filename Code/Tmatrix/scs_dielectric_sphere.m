function [scs] = scs_dieletric_sphere(k1,k2,a)
% Scattering Cross Section of a Dieleitric Sphere
%
% k1:   background wavenumber
% k2:   sphere wavenumber   
% a:    sphere radius
%
% scs:  scattering cross section of a dieletric sphere
%
% Dependencies: sbesselj, sbesselh, sbesseljp2, sbesselhp2

x1 = k1*a;
x2 = k2*a;
l = 1;
err = 1;
old = 0;
scs = 0;
while err > eps;
    j1 = sbesselj(l,x1);
    j2 = sbesselj(l,x2);
    j1p = sbesseljp2(l,x1);
    j2p = sbesseljp2(l,x2);
    h1 = sbesselh(l,x1);
    h1p = sbesselhp2(l,x1);
    num1 = j2.*j1p - j1.*j2p;
    den1 = h1.*j2p - j2.*h1p;
    num2 = x2.^2.*j2.*j1p - x1.^2.*j1.*j2p;
    den2 = x1.^2.*h1.*j2p - x2.^2.*j2.*h1p;
    t1 = abs(num1./den1).^2;
    t2 = abs(num2./den2).^2;
    tmp = (2*l+1).*(t1 + t2);
    scs = old + tmp;
    err = abs(tmp)^2/abs(scs)^2;
    old = scs;
    l = l + 1;
end
scs = 2*pi/(abs(k1)^2)*scs;