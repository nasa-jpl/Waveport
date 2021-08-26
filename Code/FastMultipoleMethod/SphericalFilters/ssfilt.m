function [f2] = ssfilt(f,L,muj,wj,K,muk,wk)
% Scalar spherical filter (interpolation or filtering)
% Assumes L <= K
%
% L:    maximum degree L<=L'
% muj:  Quadrature nodes for L'+1 points
% wj:   Quadrature weights for L'+1 ponts
% K:    maximum degree K<=K'
% muk:  Quadrature nodes for K'+1 points
% wk:   Quadrature weights for K'+1 points

% Interpolation
% f:    input scalar field f(theta_j,phi_j):    (2L'+1)x(L'+1)
% f2:   output scalar field f(theta_k,phi_k): (2K'+1)x(K'+1)
%
% Filter
% f:    input scalar field f(theta_k,phi_k):  (2K'+1)x(K'+1)
% f2:   output scalar field f(theta_j,phi_j):   (2L'+1)x(L'+1)
%
% Dependencies: sst, isst

if L > K
    error('error, L > K')
end
if L == K
    f2 = f;
    return
end
[M N] = size(f);
J = length(muj);
I = 2*(J-1) + 1;
Q = length(muk);
P = 2*(Q-1) + 1;
if M == I % interpolation mode
    if N ~= J
        error('bad array size')
    end
    filt = 0;
elseif M == P; % filter mode
    if N ~= Q
        error('bad array size')
    end
    filt = 1;
else
    error('bad array size')
end
tot = L^2 + 2*L + 1;
if ~filt % interpolation
    flm = sst(f,L,muj,wj);
    tot2 = K^2 + 2*K + 1;
    flm2 = zeros(tot2,1);
    flm2(1:tot) = flm;
    f2 = isst(flm2,K,muk);
else % filter
    flm = sst(f,K,muk,wk);
    flm2 = flm(1:tot);
    f2 = isst(flm2,L,muj);
end    