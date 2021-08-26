function [cGxy cGyx cGxz cGzx cGyz cGzy] = curlDyadicGreens(k,X,Y,Z)
% Curl of the electromagnetic free space dyadic Green's function
% Cartesian coordinates and vectors

% k:        background wavenumber
% X:        (x-x') argument
% Y:        (y-y') argument
% Z:        (z-z') argument
%           X, Y, and Z are the same size		 

% [cGxx cGyy cGzz cGxy cGxz cGyz]:		
%           6 unique dyadic components of curl G
%           Components are the same size of X, Y, and Z
%           (curl G)_xx = (curl G)_yy = (curl G)_zz = 0			
		 
% curl G =  [0	cGxy	cGxz] 
%           [cGyx	0	cGyz]		
%           [cGzx	cGzy	0]
%
%       =   [0	-d/dz	d/dy]
%           [d/dz	0	-d/dx] * exp(ik r) / (4 pi r)
%           [-d/dy	d/dx	0]

% Straight evaluation, no provisions for r=0

r = sqrt(X.^2 + Y.^2 + Z.^2);
f = (1i*k*r-1).*exp(1i*k*r)./(4*pi*(r.^3));
cGxy = -Z.*f;
cGxz = Y.*f;
cGyz = -X.*f;

% anti-symmetry
cGyx = -cGxy;
cGzx = -cGxz;
cGzy = -cGyz;

