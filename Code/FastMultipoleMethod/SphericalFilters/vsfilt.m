function [Fth2 Fphi2] = vsfilt(Fth,Fphi,L,muj,wj,K,muk,wk)
% Vector spherical filter (interpolation or filtering)
% Assumes L <= K
%
% L:    Maximum degree L<=L' to interpolate up from, or filter down to
% muj:  Quadrature nodes for L'+1 points
% wj:   Quadrature weights for L'+1 ponts
% K:    Maximum degree K<=K' to interpolate up to, or filter down from
% muk:  Quadrature nodes for K'+1 points
% wk:   Quadrature weights for K'+1 points
%
% Interpolation
% Fth, Fphi:     input fields f(theta_j,phi_j):  (2L'+1)x(L'+1)
% Fth2, Fphi2:   output fields f(theta_k,phi_k): (2K'+1)x(K'+1)
%
% Filter
% Fth, Fphi:     input scalar field f(theta_k,phi_k):  (2K'+1)x(K'+1)
% Fth2, Fphi2:   output scalar field f(theta_j,phi_j): (2L'+1)x(L'+1)
%
% Dependences: vst, ivst

if L > K
    error('error, L > K')
end
if L == K
    Fth2 = Fth;
    Fphi2 = Fphi;
    return
end
[M N] = size(Fth);
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
tot = L^2 + 2*L;
if ~filt % interpolation
    [blm clm] = vst(Fth,Fphi,L,muj,wj);
    tot2 = K^2 + 2*K;
    blm2 = zeros(tot2,1);
    clm2 = zeros(tot2,1);
    blm2(1:tot) = blm;
    clm2(1:tot) = clm;
    [Fth2 Fphi2] = ivst(blm2,clm2,K,muk);
else % filter
    [blm clm] = vst(Fth,Fphi,K,muk,wk);
    blm2 = blm(1:tot);
    clm2 = clm(1:tot);
    [Fth2 Fphi2] = ivst(blm2,clm2,L,muj);
end    
    
    
