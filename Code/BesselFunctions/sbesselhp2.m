function hp2 = sbesselhp2(n,x)
% Spherical Bessel function derivative, second type
%
%   (x h_n(x))' = h_n(x) + x h_n'(x) 
%
% n:    order of bessel function
% x:    independent variable
%
% hp2:  spherical Hankel function derivative variant
%
% Dependencies: sbesselh

hp2 = (1+n).*sbesselh(n,x) - x.*sbesselh(n+1,x);