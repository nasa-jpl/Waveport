function [w] = psdnorm(k,rmsh,lc)
% Power spectral density of a Gaussian correlated signal
%   k:      spatial frequency
%   rmsh:    root-mean-square variation
%   lc:     correlation length
    w = rmsh^2*lc*exp(-(k*lc*0.5).^2)/(2*sqrt(pi));
end