function [W, Kx, Ky] = psdOcean(Nx,Ny,Lx,Ly,u10,wind_dir_deg)
% Generate power spectral density of ocean waves based on 
% Pierson - Moskowitz (PM) wave spectrum
%
% Nx, Ny:         Number of points in x and y
% Lx, Ly:         Domain length in x and y
% u10:            Wind speed 10 meters above mean sea height
% wind_dir_deg:   Direction from which wind blows, in degrees
%               (meteorological convention)
% 
% W:              power spectral density, 2D unshifted FFT format
% Kx, Ky:         2D unshifted FFT and meshgrid format
%
% Dependencies: fftfreq

g = 9.81;                       % gravitational acceleration, (m/s^2)
alpha = 0.0081;                 % Phillips constant
wind_dir_rad = pi/180*wind_dir_deg;
fp = 0.13*g/u10;              % PM peak freqeucy (Hz)
const_f = alpha*g^2*(2*pi)^-4;  % PM scale constant
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
kx = 2*pi*fftfreq(1/dx,Nx);
ky = 2*pi*fftfreq(1/dy,Ny);
[Kx, Ky] = meshgrid(kx,ky);
K = sqrt(Kx.^2 + Ky.^2);
Theta = atan2(Ky,Kx);
F = sqrt(g*K)/2/pi;     % frequency on k-space grid, from dispersion relation
vg = g./(2*2*pi*F);     % group velocity 
cv = vg./(2*pi*K);      % change of variables 
Fk = cv.*const_f.*(F.^-5).*exp(-(5/4)*(fp./F).^4); % PM spectrum, k-space
Fk(find(K == 0)) = 0;	% zero mean
phi = (32/15)*abs(cos(0.5*(Theta-wind_dir_rad+pi)).^5);  
W = Fk.*phi;