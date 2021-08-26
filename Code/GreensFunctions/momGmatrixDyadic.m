function [Gxx Gyy Gzz Gxy Gxz Gyz] = momGmatrixDyadic(dx,kb,X,Y,Z,Xp,Yp,Zp)
% Create the dyadic MoM Green's function matrix
%
% dx:       side length of the cubic voxel (assumed constant)
% kb:       background free-space wavenumber
% X,Y,Z:    (x,y,z) coordinates, same size
% Xp,Yp,Zp: (x',y',z') coordinates, same size
%
% [Gxx Gyy Gzz Gxy Gxz Gyz]:    6 unique matrix blocks of the 3D dydaic MoM Green's
%                               Block size: MxN, M = length(X(:)), N = length(Xp(:))
%                       
% Dependencies: volintGreens, dyadicGreens

% radius of volume-equivalent sphere of cubic voxel
a = dx/(3/(4*pi))^(1/3); 

% get the volume integrated voxels
[sing delta] = volintGreens(a,kb,'dyadic');

% make primed and unprimed coordinate pairs
[XX XXP] = ndgrid(X(:),Xp(:));
[YY YYP] = ndgrid(Y(:),Yp(:));
[ZZ ZZP] = ndgrid(Z(:),Zp(:));
dX = (XX - XXP);
dY = (YY - YYP);
dZ = (ZZ - ZZP);
clear XX XXP YY YYP ZZ ZZP   % ...to save some space

% evaluate 3D dyadic Green's function 
[Gxx Gyy Gzz Gxy Gxz Gyz] = dyadicGreens(kb,dX,dY,dZ);

% multiply by the volume factor
Gxx = delta*Gxx;
Gyy = delta*Gyy;
Gzz = delta*Gzz;
Gxy = delta*Gxy;
Gxz = delta*Gxz;
Gyz = delta*Gyz;

% replace any self terms in Gxx, Gyy, Gzz with the integrated singularity
R = sqrt(dX.^2 + dY.^2 + dZ.^2);
ind = find(R==0);
Gxx(ind) = sing;
Gyy(ind) = sing;
Gzz(ind) = sing;

% replace any self terms in Gxy, Gxz, Gyz with zero
Gxy(ind) = 0;
Gxz(ind) = 0;
Gyz(ind) = 0;


