function [out1 out2 out3] = cyl2cart(rho,phi,z,Arho,Aphi,Az)
% Cylindrical to Cartesian coordinate transformation.  Angles
% in radians.
%
% rho,phi,z:    Cylindrical points
% Arho,Aphi,Az: Cylindrical vector components
%
% out1,out2,out3
%       x,y,z:        Cartesian points
%   or
%       Ax,Ay,Az:     Cartesian vector components

if checkargs(nargout,nargin)
    return
end
if nargin == 3
    out1 = rho.*cos(phi);
    out2 = rho.*sin(phi);
    out3 = z;
end
if nargin == 6
    out1 = Arho.*cos(phi) - Aphi.*sin(phi);
    out2 = Arho.*sin(phi) + Aphi.*cos(phi);
    out3 = Az;
end
    