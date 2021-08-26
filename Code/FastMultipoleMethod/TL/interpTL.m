function [tlinterp] = interpTL(theta,L,p,tlsamp,dth,thetao,M)
% Interpolation function for FMM translation operator 
%
% theta:    interpolation points, angle in radians [0 pi]
% L:        truncation degree, sum_{l=0}^L
% s:        over sampling ratio, integer
% p:        number of samples on each side of interpolation point
% Inputs from TLinterpprep:
%   tlsamp:       Translation operator sample points
%   dth:          angular sample spacing, radians
%   thetao:       stepping parameter 
%   M:            number of sample points between (0,pi), not including ends
%
% tlinterp:     interpolated translation operator values, size(theta)
%
% Dependencies: interpTLprep, SN, DM

theta = theta(:);
Nt = length(theta);
tlinterp = zeros(Nt,1);
N = M-L;
for n = 1:Nt,
    mo = floor(theta(n)/dth);
    for m=(mo-p+1):(mo+p),
        arg = theta(n)-m*dth;
        tlinterp(n) = tlinterp(n) + tlsamp(m+p)*SN(N,arg,thetao)*DM(M,arg);
    end
end
tlinterp = reshape(tlinterp,size(theta));

