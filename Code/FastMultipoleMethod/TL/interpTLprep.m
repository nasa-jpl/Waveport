function [tlsamp, dth, thetao, M, theta_samp] = interpTLprep(L,k,r,s,p)
% Preparatory function for translation operator interpolation function TLinterp.  
%
% L:    truncation degree, sum_{l=0}^L
% k:    complex wavenumber
% r:    magntude of translation distance
% s:    over sampling ratio, integer
% p:    number of samples on each side fo interpolation point
%
% tlsamp:       Translation operator sample points
% dth:          angular sample spacing, radians
% thetao:       stepping parameter
% M:            number of sample points between (0,pi), not including ends
% theta_samp:   array of theta values 
%
% Dependencies: TLth

if L < 0 || s < 0 || p < 0
    disp('bad input')
    return
end
M = s*L;
dth = 2*pi/(2*M+1);
thetao = p*dth;
theta_samp = linspace(-(p-1)*dth,pi-dth/2 + dth*p,M+2*p);
tlsamp = TLth(L,k*r,cos(theta_samp));