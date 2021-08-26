function h = rough1(Nx,Lx,rmsh,lc,type)
% Generate 1D random signal
% 
% N:        N number of points
% L:        Physical length
% rmsh:     root mean square height
% lc:       correlation length
% type:     'norm' = Gaussian power spectral density
%           'exp' = Exponential power spectral density
% 
% h:        Nx1 array of heights
%
% Dependencies: fftfreq

dx = Lx/(Nx-1);
if lc < dx
    h = rmsh*randn(Nx,1);
else
    kx = 2*pi*fftfreq(1/dx,Nx);
    if strcmp(type,'norm'),
         W = rmsh^2*lc*exp(-(kx*lc*0.5).^2)/(2*sqrt(pi));
    elseif strcmp(type,'exp'),
         W = rmsh^2*lc./(pi*(1+(kx*lc).^2));
    else
        error('bad PSD type')
    end
    gam = randn(Nx,1) + 1i*randn(Nx,1);
    gam(1) = 0;
    H = gam.*sqrt(W);
    dkx = 2*pi/Lx;
    h = (Nx*sqrt(dkx))*real(ifft(H));
end



