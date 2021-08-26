function G = momGmatrix3D(dx,kb,X,Y,Z,Xp,Yp,Zp)
% Create the 3D scalar MoM Green's function matrix
%
% dx:       side length of the cubic voxel (assumed constant)
% kb:       background free-space wavenumber
% X,Y,Z:    (x,y,z) coordinates, same size
% Xp,Yp,Zp: (x',y',z') coordinates, same size
%
% G:        3D scalar MoM Green's function matrix
%           size: MxN, M = length(X(:)), N = length(Xp(:))
%
% Dependencies: volintGreens

% radius of volume-equivalent sphere of cubic voxel
a = dx/(3/(4*pi))^(1/3); 

% get the volume integrated voxels
[sing delta] = volintGreens(a,kb,'3D');

% make primed and unprimed coordinate pairs
[XX XXP] = ndgrid(X(:),Xp(:));
[YY YYP] = ndgrid(Y(:),Yp(:));
[ZZ ZZP] = ndgrid(Z(:),Zp(:));
R = sqrt((XX - XXP).^2 + (YY - YYP).^2 + (ZZ - ZZP).^2);

% evaluate 3D scalar Green's function 
G = delta*exp(1i*kb*R)./(4*pi*R);

% replace any self terms with the integrated singularity
G(R == 0) = sing;


