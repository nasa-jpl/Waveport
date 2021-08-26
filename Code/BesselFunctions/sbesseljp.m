function jp = sbesseljp(n,x)
% Spherical Bessel function derivative
%
%	j_n'(x) = -j_(n+1)(x) + n/x * j_n(x)
%
% n:    order n
% x:    independent variable
%
% jp:   spherical bessel function derivative
%
% Dependencies: sbesselj, sbesselj2

jp = -sbesselj(n+1,x) + n.*sbesselj2(n,x);