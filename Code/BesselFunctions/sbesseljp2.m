function jp2 = sbesseljp2(n,x)
% Spherical Bessel function derivative, second type
%
%   (x j_n(x))' = j_n(x) + x j_n'(x) 
%
% n:    order of bessel function
% x:    independent variable
%
% jp2:  spherical bessel function derivative variant
%
% Dependencies: sbesselj

jp2 = (1+n).*sbesselj(n,x) - x.*sbesselj(n+1,x);
