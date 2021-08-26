function [out1 out2 out3] = cart2cyl(x,y,z,Ax,Ay,Az)
% Cartesian to cylindrical coordinate transformation.  Angles
% in radians.
%
% x,y,z:        Cartesian points
% Ax,Ay,Az:     Cartesian vector components
%
% out1,out2,out3
%       rho,phi,z:    Cylindrical points
%    or 
%       Arho,Aphi,Az: Cylindrical vector components

if checkargs(nargout,nargin)
    return
end
if nargin == 3
    out1 = sqrt(x.^2 + y.^2);
    out2 = atan2(y,x);
    out3 = z;
end
if nargin == 6
    phi = atan2(y,x);
    out1 = Ax.*cos(phi) + Ay.*sin(phi);
    out2 = -Ax.*sin(phi) + Ay.*cos(phi);
    out3 = Az;
end
