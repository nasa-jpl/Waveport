function hp = sbesselhp(n,x)
% Spherical Hankel function derivative
%
%	h_n'(x) = -h_(n+1)(x) + n/x * h_n(x)
%
% n:    order n
% x:    independent variable
%
% hp:   spherical Hankel function derivative
%
% Dependencies: sbesselh

hp = -sbesselh(n+1,x) + (n./x).*sbesselh(n,x);

