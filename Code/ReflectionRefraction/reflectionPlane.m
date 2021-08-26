function [rx ry rz theta l1 l2] = reflectionPlane(x1,y1,z1,x2,y2,z2,zo);
% Reflection point on a plane
%
% x1,y1,z1      Coordinates of exterior point, r1, centered on circle
% x2,y2,z1      Coordinates of interior point, r2, centered on circle
% zo            z-level of the plane
%
% rx,ry,rz      Coordinates of refraction point
% theta         Reflection angle
% l1,l2         Lengths of the vectors from each point to the reflection point

% note bad z-coordinates
ind = or((z1 <= 0),(z2 <= 0));

% compute parameters
h1 = z1 - zo;
h2 = z2 - zo;
L = sqrt((x2-x1).^2 + (y2-y1).^2);
uhat_x = (x2-x1)./L;
uhat_y = (y2-y1)./L;

% zero out points on top of each other
indL = (L==0);
uhat_x(indL) = 0;
uhat_y(indL) = 0;

% solution to u1
u1 = h1.*L./(h2 + h1);
u2 = L - u1;

% compute reflection point
rx = x1 + u1.*uhat_x;
ry = y1 + u1.*uhat_y;
rz = zo*ones(size(x1));

% compute other quantities
theta = atan2(u1,h1);
l1 = sqrt(u1.^2 + h1.^2);
l2 = sqrt(u2.^2 + h2.^2);

% nan-out bad z-coordinates
rx(ind) = nan;
ry(ind) = nan;
rz(ind) = nan;
theta(ind) = nan;
l1(ind) = nan;
l2(ind) = nan;


