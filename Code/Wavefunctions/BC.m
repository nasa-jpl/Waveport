function [Bth, Bphi, Cth, Cphi] = BC(L,theta,phi,normstr)
% Vector spherical harmonics, defaults to partial normalization
%
% B_lm(theta,phi) = N_lm [theta_hat d/dtheta P_l^m(cos(theta)) + 
%                         phi_hat im/sin(theta) P_l^m(cos(theta))]e^(im phi)
%
% C_lm(theta,phi) = N_lm [theta_hat im/sin(theta) P_l^m(cos(theta)) - 
%                         phi_hat d/dtheta P_l^m(cos(theta))]e^(im phi)                       
%
% N_lm = (-1)^m/sqrt(l(l+1))*sqrt((2l+1)/4pi*(l-m)!/(l+m)!)
%
% L:                    Maximum harmonic degree L
% theta, phi:           Spherical angles in radians, any size
% normstr:              [optional] 'norm' for 1/sqrt(l(l+1)) full normalization
%
% Bth,Bphi,Cth,Cphi:    Vector spherical harmonics components 
%                       evaluated at (theta,phi)
%                       Dimensions: length(theta(:)) x N
%                       N = L^2 + 2*L
%
% Dependencies: ind2lm, Plmp2, mPlmsin

N = L^2 + 2*L;
theta = theta(:)';
phi = phi(:)';

% d/dtheta P_l^m
dP = Plmp2(L,theta).';

% m P_l^m/sin(theta)
mPsin = mPlmsin(L,theta).';

% remaining normalizaton and exponential
[l, m] = ind2lm((1:N)');
[Phi,M] = ndgrid(phi,m);
ex = 1./sqrt(2*pi).*exp(1i*M.*Phi);
if nargin == 4 
    if ~strcmp(normstr,'norm')
        error('bad string')
    end
    [~,L] = ndgrid(phi,l);
    ex =  ex.*1./sqrt(L.*(L+1));
end
dP = dP.*ex;
mPsin = 1i*mPsin.*ex;

% B and C components theta_hat, phi_hat
Bth = dP;
Bphi = mPsin;
Cth = mPsin;
Cphi = -dP;

