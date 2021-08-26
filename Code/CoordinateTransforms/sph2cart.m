function [out1 out2 out3] = sph2cart(r,th,phi,Ar,Ath,Aphi)
% Cartesian to spherical coordinate transformation.  Angles
% in radians.
%
% rho,phi,z:        Spherical points
% Ar,Ath,Aphi:      Spherical vector components
%
% out1,out2,out3
%       x,y,z:        Cartesian points
%    or 
%       Ax,Ay,Az:     Cartesian vector components

if checkargs(nargout,nargin)
    return
end
if nargin == 3
    out1 = r.*sin(th).*cos(phi);
    out2 = r.*sin(th).*sin(phi);
    out3 = r.*cos(th);
end
if nargin == 6
    out1 = Ar.*sin(th).*cos(phi) + Ath.*cos(th).*cos(phi) - Aphi.*sin(phi);
    out2 = Ar.*sin(th).*sin(phi) + Ath.*cos(th).*sin(phi) + Aphi.*cos(phi);
    out3 = Ar.*cos(th) - Ath.*sin(th);
end
    