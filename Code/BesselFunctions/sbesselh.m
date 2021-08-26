function h = sbesselh(n,x)
% Spherical Hankel function
%
%   h_n(x) = sqrt(pi/(2x)) H_(n+1/2)(x)
%
% n:    order n
% x:    independent variable
%
% h:    spherical Hankel function order n evaluated at x

h = sqrt(pi./(2*x)).*besselh(n+0.5,x);
