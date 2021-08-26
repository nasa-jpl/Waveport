function [ylm] = sphericalY(L,theta,phi,str)
% Fully normalized spherical harmonics
%
% Y_lm(theta,phi) = sqrt((2l+1)/4pi*(l-m)!/(l+m)!)
%                       *P_l^m(cos(theta)) e^(i m phi)
%
%                   (-1)^m included in P_l^m
%
% L:            Maximum harmonic degree L
% theta, phi:   Spherical angles in radians
% str:          [optional] 'mono' to include monopole term
%
% ylm:          Spherical harmonics evaluated at (theta,phi)
%               Dimensions: length(theta(:)) x N
%               N = L^2 + 2*L (N = L^2 + 2*L + 1 with monopole)
%
% Dependencies: ind2lm, Plm

N = L^2 + 2*L + 1;
[l, m] = ind2lm((1:N),'mono');
[Phi,M] = ndgrid(phi(:),m);
ylm = (1/sqrt(2*pi))*Plm(L,cos(theta(:))).'.*exp(1i*M.*Phi);
if nargin == 4
    if ~strcmp(str,'mono')
        error('bad string')
    end
    return
else
    ylm = ylm(:,2:end);
end

