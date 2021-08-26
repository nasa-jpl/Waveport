function [flm] = sst(f,L,muj,wj)
% Forward scalar spherical transform
%
% f:    [IxJ] sampled spherical function
%       I = 2*L'+1, J = L'+1
%       phi_i = 2*pi*(0:(I-1))/I
%       theta = arccos(mu_j)
% L:    maximum degree L for output coefficients flm (L' >= L)
% muj:  J Gaussian quadrature nodes
% wj:   J Gaussian quadrature weights
%
% flm:  [Nx1], N = L^2+2*L+1, spectral coefficients, linearly indexed
%
% Dependencies: Plm

[I, J] = size(f);
tot = L^2 + 2*L + 1;
flm = zeros(tot,1);
% phi transform
fmth = (sqrt(2*pi)/I)*fft(f,[],1);
% Legendre polynomials evaluatd at muj
plm = Plm(L,muj);
% theta quadrature integral
wj = wj(:).';
for l=0:L,
for m=-l:l,
    ind = lm2ind(l,m,'mono');
    if m >= 0
        ind2 = m + 1;
    else
        ind2 = I + m + 1;
    end
    flm(ind) = sum(fmth(ind2,:).*plm(ind,:).*wj); 
end
end

