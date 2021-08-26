function [w] = psdexp(k,rmsh,lc)
% Power spectral density of exponential correlated signal
%   k:      spatial frequency
%   rmsh:    root-mean-square variation
%   lc:     correlation length
    w = rmsh^2*lc./(pi*(1+(k*lc).^2));
end