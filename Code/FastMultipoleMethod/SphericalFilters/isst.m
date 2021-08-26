function [f] = isst(flm,L,muj)
% Inverse scalar spherical transform
%
% flm:  tot = L^2+2*L+1 spectral coefficients
% L:    maximum degree L
% muj:  J Gaussian quadrature nodes, muj = cos(theta_j)
%       
% f:    [IxJ] gridded spherical function
%       I = 2*L'+1, J = L'+1, where J = length(muj)
%
% Dependencies: Plm

J = length(muj);
I = 2*(J-1) + 1;
% Legendre polynomials evaluatd at muj
plm = Plm(L,muj);
% inverse Legendre transform
fmth = zeros(I,J);
for m = -L:L,
    if m >= 0
        ind2 = m + 1;
    else
        ind2 = I + m + 1;
    end    
    for l = abs(m):L,
        ind = lm2ind(l,m,'mono');
        fmth(ind2,:) = fmth(ind2,:) + flm(ind)*plm(ind,:);
    end            
end
% inverse Fouriner transform in phi
f = (I/sqrt(2*pi))*ifft(fmth,[],1);


