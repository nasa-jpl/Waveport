function [Gxx Gyy Gzz Gxy Gxz Gyz] = dyadicGreens(k,X,Y,Z)
% Electromagnetic free space dyadic Green's function in 
% Cartesian coordinates. No provisions for r=0.
%
% k:	background wavenumber
% X:	(x-x') 
% Y:	(y-y') 
% Z:	(z-z') 
%       X, Y, and Z are the same size		 
%
% [Gxx Gyy Gzz Gxy Gxz Gyz]: Unique matrix entries of G(r,r')

r = sqrt(X.^2 + Y.^2 + Z.^2);
c1 = -1./(r.^2)-3*1i./k./(r.^3)+3./(k.^2)./(r.^4);
c2 = 1i./k./r-1./(k.^2)./(r.^2);
c3 = exp(1i*k.*r)./(4*pi*r);
Gxx = (1+c1.*X.^2+c2).*c3;
Gyy = (1+c1.*Y.^2+c2).*c3;
Gzz = (1+c1.*Z.^2+c2).*c3;
Gxy = c1.*X.*Y.*c3;
Gxz = c1.*X.*Z.*c3;
Gyz = c1.*Y.*Z.*c3;


