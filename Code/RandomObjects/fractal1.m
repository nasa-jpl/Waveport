function z = fractal1(Nx,Lx,H,sig_o,x_o,type)
% 1D Fractal Brownian profile
%
% Nx:       Number of points in the profile
% Lx:       Profile length
% H:        Hurst exponent between [0,1];
% sig_o:    RMS height \sigma_o at scale x_o (or RMS deviation v_o at lag \Delta x_o)
% x_o:      Length scale x_o at which sig_o applies (or lag \Delta x_o at which v_o applies)
% type:     'rms' for RMS height (sig_o, x_o) 
%           'dev' for RMS deviation (v_o, \Delta x_o)
%
% z:        1D height profile
%
% Notes: Based on fbm1d.m (Kroese, et. al, 2015)    
    
if H < 0 || H > 1
    error('Bad Hurst exponent')
end
if strcmp(type,'rms')
    typein = 1;
elseif strcmp(type,'dev')
    typein = 0;
else
    error('Bad type')
end
% covariance
R = nan(Nx+1,1);
R(1) = 1;
ind = 1:Nx;
R(ind+1) = 0.5*((ind+1).^(2*H) - 2*ind.^(2*H) + (ind-1).^(2*H));
R = [R; R(end-1:-1:2)]; % load the circulant matrix
lambda = real(fft(R))/(2*Nx); % eigenvalues
gam = randn(2*Nx,1) + 1i*randn(2*Nx,1); % random seed
z = real(fft(gam.*sqrt(lambda)));
z = Nx^(-H)*cumsum(z(1:Nx+1)); 
z = z(1:Nx);
z = z-mean(z); % zero mean
const = sig_o*(Lx/x_o)^H;
% normalize according to RMS height or RMS variation
if typein
    z = const*z/std(z);
else
    z = const*z;
end

    
    