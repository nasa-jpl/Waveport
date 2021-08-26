function At = evolveOcean(W,Kx,Ky,gam,t,type)
% Evolve an ocean surface at time t given the spectrum and time
% Creamer2 nonlinear surface is an option.  
%
% W:            Ocean spectrum from psdOcean
% gam:          One realization of 2D complex standard normal
% Kx,Ky:        Spatial frequency matrices from psdOcean
% t:            Time (s)
% type:         'creamer2' for nonlinear surface
% 
% At:           Ocean surface at time t
%
% Dependencies: psdOcean

g = 9.81;           % gravitational acceleration, (m/s^2)
K = sqrt(Kx.^2 + Ky.^2);
ww = sqrt(g*K);     % dispersion relation
dkx = Kx(1,2) - Kx(1,1);
dky = Ky(2,1) - Ky(1,1);
Atk = gam.*sqrt(dkx*dky*W).*exp(-1i*ww*t); % linear complex amplitudes
[Nx Ny] = size(W);
At = (Nx*Ny)*real(ifft2(Atk)); % linear surface height
if nargin == 6 && strcmp(type,'creamer2')
   % Hilbert Transform
   htk_x = -1i*(Kx./K).*Atk;
   htk_y = -1i*(Ky./K).*Atk;
   indInf = find(K == 0);
   htk_x(indInf) = 0;
   htk_y(indInf) = 0;
   htx_x = real(ifft2(htk_x));
   htx_y = real(ifft2(htk_y));
   % Creamer 2 spectrum
   Ct2k = -(Kx.^2)./(2*K).*fft2(htx_x.^2) ...
         -(Kx.*Ky)./(K).*fft2(htx_x.*htx_y) ...
         -(Ky.^2)./(2*K).*fft2(htx_y.^2);
   Ct2k(indInf) = 0;  
   % Radial lowpass filter
   kcutoff = 5;
   mu = 8;
   Bk = 1./(1+(K/kcutoff).^(mu));
   At = Nx*Ny*real(ifft2(Atk + Bk.*Ct2k));
end
   