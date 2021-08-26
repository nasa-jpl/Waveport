function [h] = rough3(Nx,Ny,Nz,Lx,Ly,Lz,rmsh,lcx,lcy,lcz,type)
% Generate 3D random volume
%
%   Nx,Ny,Nz:       number of points in x, y, and z
%   Lx,Ly,Lz:       physical extent of x, y, and z
%   rmsh:           root mean square height
%   lcx,lcy,lcz:    correlation length
%   type:           'norm' = Gaussian power spectral density
%                   'exp' = Exponential power spectral density
%
%   h:              [Ny,Nx,Nz] meshgrid array of values
%
% Dependencies: fftfreq

if strcmp(type,'norm')
    typesw = 1;
elseif strcmp(type,'exp')
    typesw = 0;
else
   error('bad PSD type')
end    
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
dz = Lz/(Nz-1);
dd = [dy dx dz];
Ns = [Ny Nx Nz];
Ls = [Ly Lx Lz];
Lcs = [lcy lcx lcz];
undersamp = [(lcy<dy) (lcx<dx) (lcz<dz)];
% all correlation lengths are well sampled
if sum(undersamp) == 0
    kx = 2*pi*fftfreq(1/dx,Nx);
    ky = 2*pi*fftfreq(1/dy,Ny);
    kz = 2*pi*fftfreq(1/dz,Nz);
    [Kx Ky Kz] = meshgrid(kx,ky,kz);
    if typesw,
        W = (rmsh^2*lcx*lcy*lcz/8/(pi^(3/2)))*exp(-lcx^2*Kx.^2/4-lcy^2*Ky.^2/4-lcz^2*Kz.^2/4);
    else
        W = (rmsh^2*lcx*lcy*lcz)./(3*pi*(1+(Kx*lcx).^2+(Ky*lcy).^2+(Kz*lcz).^2).^2);
    end
    gam = randn(Ny,Nx,Nz) + 1i*randn(Ny,Nx,Nz);
    gam(1,1,1) = 0;
    W = gam.*sqrt(W);
    dkx = 2*pi/Lx;
    dky = 2*pi/Ly;
    dkz = 2*pi/Lz;
    h = (Nx*Ny*Nz*sqrt(dkx*dky*dkz))*real(ifftn(W));
% one correlation length is less than its sample rate
% shiftdim to put this in the last dimension
% use 2D PSD in dimensions 1 and 2, then shift back
elseif sum(undersamp) == 1
    dim = find(undersamp);
    shift = -dim;
    Ns = circshift(Ns,shift);
    dd = circshift(dd,shift);
    Ls = circshift(Ls,shift);
    Lcs = circshift(Lcs,shift);
    kx = 2*pi*fftfreq(1/dd(1),Ns(1));
    ky = 2*pi*fftfreq(1/dd(2),Ns(2));
    [Kx Ky] = meshgrid(kx,ky);
    Kx = repmat(Kx,1,1,Ns(3));
	Ky = repmat(Ky,1,1,Ns(3));
    lcx = Lcs(1); 
    lcy = Lcs(2);
    if typesw
        W = (rmsh^2*lcx*lcy/4/pi)*exp(-lcx^2*Kx.^2/4-lcy^2*Ky.^2/4);
    else
        W = (rmsh^2*lcx*lcy)./(2*pi*(1+(Kx*lcx).^2+(Ky*lcy).^2).^(3/2));
    end
    gam = randn(Ns) + 1i*randn(Ns);
    gam(1,1,:) = 0;
    W = gam.*sqrt(W);
    dkx = 2*pi/Ls(1);
    dky = 2*pi/Ls(2);
    h = (Nx*Ny*sqrt(dkx*dky))*real(ifft(ifft(W,[],1),[],2));  
    h = shiftdim(h,3+shift);
% two correlation lengths are less than their samples rates
% shiftdim to put one sampled dimension in dimension 1
% use 1D PSD in dimension 1, then shift back
elseif sum(undersamp) == 2
    dim = find(~undersamp);
    shift = 1-dim;
    Ns = circshift(Ns,shift);
    dd = circshift(dd,shift);
    Ls = circshift(Ls,shift);
    Lcs = circshift(Lcs,shift);
    kx = 2*pi*fftfreq(1/dd(1),Ns(1));
    K = repmat(kx,1,Ns(2),Ns(3));
    lc = Lcs(1);
    if typesw
        W = (rmsh^2*lc/2/sqrt(pi))*exp(-lc^2*K.^2/4);
    else 
        W = (rmsh^2*lc)./(pi*(1+(K*lc).^2));
    end
    gam = randn(Ns) + 1i*randn(Ns);
    gam(1,:,:) = 0;
    W = gam.*sqrt(W);
    dk = 2*pi/Ls(1);
    h = (Ns(1)*sqrt(dk))*real(ifft(W,[],1));  
    h = shiftdim(h,3+shift);
% all correlation lengths less than sample rate
else 
    h = rmsh*randn(Ny,Nx,Nz);
end