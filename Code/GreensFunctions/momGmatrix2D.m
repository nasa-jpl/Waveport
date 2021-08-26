function G = momGmatrix2D(dx,kb,X,Y,Xp,Yp)
% Create the 2D scalar MoM Green's function matrix
%
% dx:       side length of the square voxel (assumed constant)
% kb:       background free-space wavenumber
% X,Y:      (x,y) coordinates, same size
% Xp,Yp:    (x',y') coordinates, same size
%
% G:        2D scalar MoM Green's function matrix
%           size: MxN, M = length(X(:)), N = length(Xp(:))
%
% Dependencies: volintGreens

% radius of area-equivalent circle of square voxel
a = dx/sqrt(pi);

% get the volume integrated voxels
[sing delta] = volintGreens(a,kb,'2D');

% make primed and unprimed coordinate pairs
[XX XXP] = ndgrid(X(:),Xp(:));
[YY YYP] = ndgrid(Y(:),Yp(:));
Rho = sqrt((XX - XXP).^2 + (YY - YYP).^2);

% evaluate 2D scalar Green's function 
G = (delta*1i/4)*besselh(0,kb*Rho);

% replace any self terms with the integrated singularity
G(Rho == 0) = sing;