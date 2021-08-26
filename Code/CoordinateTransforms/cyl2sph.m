function [out1 out2 out3] = cyl2sph(rho,phi,z,Arho,Aphi,Az)
% Cyindrical to spherical coordinate transformation.  Angles
% in radians.
%
% rho,phi,z:        Cyindrical points
% Arho,Aphi,Az:     Cyindrical vector components
%
% out1,out2,out3
%       r,th,phi:       Spherical points
%    or 
%       Ar,Ath,Aphi:	Spherical vector components

if checkargs(nargout,nargin)
    return
end
if nargin == 3
    out1 = sqrt(rho.^2 + z.^2);
    out2 = atan2(rho,z);
    out3 = phi;
end
if nargin == 6
    th = atan2(rho,z);
    out1 = Arho.*sin(th) + Az.*cos(th);
    out2 = Arho.*cos(th) - Az.*sin(th);
    out3 = Aphi;
end
