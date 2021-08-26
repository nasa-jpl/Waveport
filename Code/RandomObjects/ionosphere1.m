function [Ne V kx Lo Lb v1 v2] = ionosphere1(Nx,Lx,ave,pctrms,varargin)
% 1D profile of ionosphere irregularity using 2-parameter model
% 
% N:            N number of points
% L:            Profile length (m)
% ave:          Average value of electron density or TEC
% pctrms:       Percent (%) RMS of electron density or TEC deviation
% varargin:     [optional] Lo, Lb, v1, v2
%               Lo, Lb: Outerscale and breakscale (m)
%               v1, v2: log-slope parameters
% 
% Ne:           Nx1 array of ionosphere irregularity profile (electron density, Ne, or TEC)
% V,kx:         [optional] Spectral density and wavenumber 
% Lo,Lb,v1,v2:  [optional] will output default values if needed
%
% Dependencies: fftfreq

dx = Lx/(Nx-1);
kx = 2*pi*fftfreq(1/dx,Nx);

% spectral length scales and log-slope parameters
% defaults
Lo = 10e3;      % outerscale, 10 km
Lb = 500;       % breakscale, 500 m
v1 = 3.5/2;     % log-slope parameter 1
v2 = 5.5/2;     % log-slope parameter 2

% user defined parameters
if length(varargin)>=1 && ~isempty(varargin{1}), Lo = varargin{1}; end
if length(varargin)>=2 && ~isempty(varargin{2}), Lb = varargin{2}; end
if length(varargin)>=3 && ~isempty(varargin{3}), v1 = varargin{3}; end
if length(varargin)==4 && ~isempty(varargin{4}), v2 = varargin{4}; end

% wavenumbers
ko = 2*pi/Lo;   % outerscale wavenumber
kb = 2*pi/Lb;   % breakscale wavenumber

% initialize spectral density, V
V = zeros(Nx,1);

% compute V at lower wavenumbers
ind = (abs(kx) <= kb);
V(ind) = (pi/(v1-1))./((kx(ind).^2+ko^2).^(v1-1)) ...
    - pi/((kb^2+ko^2)^(v1-1))*(1/(v1-1)-1/(v2-1));

% compute V at upper wavenumbers
ind = (abs(kx) > kb);
V(ind) = (pi/(v2-1))*((ko.^2+kb^2)^(v2-v1))./((kx(ind).^2+ko^2).^(v2-1));

% create unnormalized profile
gam = randn(Nx,1) + 1i*randn(Nx,1);
gam(1) = 0;
H = gam.*sqrt(V);
dkx = 2*pi/Lx;
Ne = (Nx*sqrt(dkx))*real(ifft(H));

% normalize by RMS and apply scale factors
Ne = Ne/rms(Ne);
Ne = ave*(pctrms/100*Ne + 1);