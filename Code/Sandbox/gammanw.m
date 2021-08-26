function [g] = gammanw(sz)
% Frequency domain complex standard normal distribution
%
% sz:   Output dimensions [N1, N2, ...]
%
% g:    Complex standard normal with DFT conjugate symmetry

g = (1/sqrt(prod(sz)))*fftn(randn(sz));
g(1) = 0;  % zero the DC term
