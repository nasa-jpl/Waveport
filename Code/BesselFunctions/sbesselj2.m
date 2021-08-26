function j = sbesselj2(n,x)
% Spherical Bessel function variant
%
%   j_n(x)/x
%
% n:    order n
% x:    independent variable
%
% j:    spherical bessel function order n evaluated at x, divided by x

j = sbesselj(n,x)./x;
thresh = 1e-5;
j(and((x<thresh),(n==1))) = 1/3 - x(and((x<thresh),(n==1))).^2/30;
j(and((x==0),(n>1))) = 0;