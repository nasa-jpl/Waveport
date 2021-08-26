function j = sbesselj(n,x)
% Spherical Bessel function
%
%   j_n(x) = sqrt(pi/(2x)) J_(n+1/2)(x)
%
% n:    order n
% x:    independent variable
%
% j:    spherical bessel function order n evaluated at x

j = sqrt(pi./(2*x)).*besselj(n+0.5,x);
thresh = 1e-5;
j(and((x<thresh),(n==0))) = 1 - x(and((x<thresh),(n==0))).^2/6;
j(and((x==0),(n~=0))) = 0;


