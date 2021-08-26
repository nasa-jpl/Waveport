function [th, phi] = discoball(nlat)
% Approximate uniform distribution of points on a sphere
% using a disco ball method.
%
% nlat:         number of latitude lines
%
% th, phi:      spherical angles of points, in radians.     

th = linspace(0,pi,nlat)';
dth = th(2)-th(1);
xx = [];
for n=1:nlat,
    rad = sin(th(n));    
    nlon = round(2*pi*rad/dth)+1;
    dlon = 2*pi/nlon;
    lon = linspace(0,2*pi-dlon,nlon)';
    xx = [xx; th(n)*ones(nlon,1) lon];
end
th = xx(:,1);
phi = xx(:,2);

