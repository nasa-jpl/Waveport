function [out1 out2 out3] = sph2cyl(r,th,phi,Ar,Ath,Aphi)
% Spherical to cylinder coordinate transformation.  Angles
% in radians.
%
% r,th,phi:        Spherical points
% Ar,Ath,Aphi:      Spherical vector components
%
% out1,out2,out3
%       rho,phi,z:    	Cylinder points
%    or 
%       Arho,Aphi,Az: 	Cylinder vector components

if checkargs(nargout,nargin)
    return
end
if nargin == 3
    out1 = r.*sin(th);
    out2 = phi;
    out3 = r.*cos(th);
end
if nargin == 6
    out1 = Ar.*sin(th) + Ath.*cos(th);
    out2 = Aphi;
    out3 = Ar.*cos(th) - Ath.*sin(th);
end
