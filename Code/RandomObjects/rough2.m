function [h] = rough2(Nx,Ny,Lx,Ly,rmsh,lcx,lcy,type)
% Generate 2D random rough surface 
%
%  Nx,Ny:       number of points in x and y
%  Lx,Ly:       physical extent of x and y
%  rmsh:        root mean square height
%  lcx,lcy:     correlation lengths in x and y
%  type:        'norm' = Gaussian power spectral density
%               'exp' = Exponential power spectral density
%
%  h:           [Ny,Nx] meshgrid array of heights
%
%  Dependencies: fftfreq

if strcmp(type,'norm')
    typesw = 1;
elseif strcmp(type,'exp')
    typesw = 0;
else
   error('bad PSD type')
end    
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
dd = [dx dy];
Ns = [Ny Nx];
Ls = [Ly Lx];
Lcs = [lcy lcx];
undersamp = [(lcy<dy) (lcx<dx)];
% all correlation lengths are well sampled
if sum(undersamp) == 0
    kx = 2*pi*fftfreq(1/dx,Nx);
    ky = 2*pi*fftfreq(1/dy,Ny);
    [Kx Ky] = meshgrid(kx,ky);
    if typesw
        W = (rmsh^2*lcx*lcy/4/pi)*exp(-lcx^2*Kx.^2/4-lcy^2*Ky.^2/4);
    else
        W = (rmsh^2*lcx*lcy)./(2*pi*(1+(Kx*lcx).^2+(Ky*lcy).^2).^(3/2));
    end
    gam = randn(Ns) + 1i*randn(Ns);
    W = gam.*sqrt(W);
    dkx = 2*pi/Lx;
    dky = 2*pi/Ly;
    h = Nx*Ny*sqrt(dkx*dky)*real(ifft2(W));  
% one correlation length is less than the sample rate
elseif sum(undersamp) == 1
    shift = find(~undersamp)-1;
    Ns = circshift(Ns,shift);
    dd = circshift(dd,shift);
    Ls = circshift(Ls,shift);
    Lcs = circshift(Lcs,shift);
    kx = 2*pi*fftfreq(1/dd(1),Ns(1));
    K = repmat(kx,1,Ns(2));
    lc = Lcs(1);
    if typesw
        W = (rmsh^2*lc/2/sqrt(pi))*exp(-lc^2*K.^2/4);
    else 
        W = (rmsh^2*lc)./(pi*(1+(K*lc).^2));
    end
	gam = randn(Ns) + 1i*randn(Ns);
    gam(1,:) = 0;
    W = gam.*sqrt(W);
    dk = 2*pi/Ls(1);
    h = (Ns(1)*sqrt(dk))*real(ifft(W,[],1));  
    h = shiftdim(h,shift);
% both correlation lengths less than their sample rates
else
    h = rmsh*randn(Ny,Nx);
end


           
