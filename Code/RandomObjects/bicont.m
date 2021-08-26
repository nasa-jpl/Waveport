function [Theta S] = bicont(X,Y,Z,lave,b,fv,N)
% Bicontinuous random media 
%
% X,Y,Z     [any size] Cartesian sample points
% lave      Average length scale of heterogeneity
% b         Shape parameter
% fv        Volume fraction between 0 and 1
% N         [optional] Number of plane waves
%
% Theta     Bicontinuous indicator function
% S         [optional] Fluctuating field
%
% Dependencies: randsphere,sph2cart

if nargin == 6
    N = 10^3;
end

% scale and shape parameters
kave = 2*pi/lave;
shape = b+1;
scale = kave/(b+1);
alp = erfinv(1-2*fv);

% compute random wave vectors
k_r = gamrnd(shape,scale,N,1);
[k_th k_ph] = randsphere(N);
[kx ky kz] = sph2cart(k_r,k_th,k_ph);
phi = 2*pi*rand(N,1);

% create fluctuating field
S = zeros(size(X));
for n=1:N
    S = S + cos(kx(n)*X + ky(n)*Y + kz(n)*Z + phi(n));
end
S = S/sqrt(N);

% indicator function from level cut
Theta = zeros(size(S));
Theta(S > alp) = 1;


