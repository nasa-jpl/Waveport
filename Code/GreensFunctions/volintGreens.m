function [sing delta] = volintGreens(a,kb,str)
% Volume-intergrated scalar and dyadic Green's function
%
% a:            radius of circular or spherical voxel
% kb:           background wavenumber
% str:          '2D'        : 2D scalar VIE
%               '3D'        : 3D scalar VIE
%               'dyadic'    : 3D vector VIE
%
% sing:         Green's function integrated at the singular point
% delta:        Discrete surface or volume element for non-singular points

% singular point
if strcmp(str,'2D') % 2D scalar VIE
    sing = -1./kb.^2 + (1i*pi/2)*(a./kb).*besselh(1,kb.*a);
elseif strcmp(str,'3D') % 3D scalar VIE
    sing = (1./kb.^2).*(-1 + (1 - 1i*kb.*a).*exp(1i*kb.*a));
elseif strcmp(str,'dyadic') % 3D vector VIE
    sing = (1./kb.^2).*(-1 + (2/3)*(1 - 1i*kb.*a).*exp(1i*kb.*a));
else
    error('Bad string')
end

% Discrete surface or volume element for 2D or 3D 
if strcmp(str,'2D')
    delta = 2*pi*(a./kb).*besselj(1,kb.*a);
else
    delta = 4*pi*(a./kb.^2).*(sin(kb.*a)./(kb.*a) - cos(kb.*a));
end




