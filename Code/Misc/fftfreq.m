function [f] = fftfreq(fs,N)
% Two-sided FFT frequencies
%
% fs:   sample rate
% N:    number of samples
%
% f:    fft frequencies (same units as fs)

if mod(N,2) == 0    % N even
    seg = 1:(N/2-1);
    f = [0 seg -N/2 -fliplr(seg)]'*fs/N;    
else                % N odd
    seg = 1:((N-1)/2);
    f = [0 seg -fliplr(seg)]'*fs/N;
end
