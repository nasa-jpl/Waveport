function [th, phi] = randsphere(N)
% Uniformly random points on a sphere
%
% N:        number of points
%
% th, phi:  spherical angles of points, in radians.     

th = acos(2*rand(N,1) - 1);
phi = 2*pi*rand(N,1);

