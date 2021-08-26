function [out1 out2 out3] = cart2sph(x,y,z,Ax,Ay,Az)
% Cartesian to spherical coordinate transformation.  Angles
% in radians.
%
% x,y,z:        Cartesian points
% Ax,Ay,Az:     Cartesian vector components
%
% out1,out2,out3
%       r,theta,phi:   Spherical points
%    or 
%       Ar,Ath,Aphi: Spherical vector components

if checkargs(nargout,nargin)
    return
end
if nargin == 3
    out1 = sqrt(x.^2 + y.^2 + z.^2);
    out2 = atan2(sqrt(x.^2 + y.^2),z);
    out3 = atan2(y,x);
end
if nargin == 6
    th = atan2(sqrt(x.^2 + y.^2),z);
    phi = atan2(y,x);
    out1 = Ax.*sin(th).*cos(phi) + Ay.*sin(th).*sin(phi) + Az.*cos(th);
    out2 = Ax.*cos(th).*cos(phi) + Ay.*cos(th).*sin(phi) - Az.*sin(th);
    out3 = -Ax.*sin(phi) + Ay.*cos(phi);
end
