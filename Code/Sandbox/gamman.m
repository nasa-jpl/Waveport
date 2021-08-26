function [g] = gamman(sz)
% Complex standard normal distribution
%
% sz:   Output dimensions [N1, N2, ...]
%
% g:    Complex standard normal

g = (1/sqrt(2))*(randn(sz) + 1i*randn(sz));

