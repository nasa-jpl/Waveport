function [Fth2 Fphi2] = fvsfilt1(Fth,Fphi,L,K,Am,Bm)
% Fast vector spherical filter (interpolation or filtering)
% Based on precomputed matrix-vector multipy.
% Assumes L <= K. Auto-detects filter or interplation based on size of Fth,
% Fphi. Am, Bm are computed with fvsfilt1AmBm for filter or interpolation.
%
% L:        Lower degree harmonic (to interpolate up from, or filter down to) 
% K:        Upper degree harmonic (to interpolate up to, or filter down from)
% Am,Bm:    [JxQxI] matrix when interpolating from L to K
%           [QxJxP] matrix when filtering from K to L
%
% Interpolation
% Fth, Fphi:     input fields f(theta_j,phi_j):  (2L+1)x(L+1)
% Fth2, Fphi2:   output fields f(theta_k,phi_k): (2K+1)x(K+1)
%
% Filter
% Fth, Fphi:     input scalar field f(theta_k,phi_k):  (2K+1)x(K+1)
% Fth2, Fphi2:   output scalar field f(theta_j,phi_j): (2L+1)x(L+1)
%
% Dependences: fvsfilt1AmBm

if L > K, error('error, L > K'), end
if L == K
    Fth2 = Fth;
    Fphi2 = Fphi;
    return
end
[M N] = size(Fth);
I = 2*L+1;
J = L+1;
P = 2*K+1;
Q = K+1;
if M == I % interpolation 
    if N ~= J, error('bad array size'), end
    filt = 0;
    Fth2 = zeros(P,Q);
    Fphi2 = zeros(P,Q);
    const = P/I;
elseif M == P; % filter 
    if N ~= Q, error('bad array size'), end
    filt = 1;
    Fth2 = zeros(I,J);
    Fphi2 = zeros(I,J);
    const = I/P;
else
    error('bad array size')
end

% compute fft in phi
fthm = const*fft(Fth,[],1);
fphim = const*fft(Fphi,[],1);
% loop over m at L
for m = -L:L,
    if m >= 0
        ind1 = m + 1;
        ind2 = m + 1;
    else
        ind1 = I + m + 1;
        ind2 = P + m + 1;
    end
    % perform matrix vector multiply for each m
    if filt % filtering
        Fth2(ind1,:)  = fthm(ind2,:)*Am(:,:,ind1) + fphim(ind2,:)*Bm(:,:,ind1);
        Fphi2(ind1,:) = -fthm(ind2,:)*Bm(:,:,ind1) + fphim(ind2,:)*Am(:,:,ind1);      
    else % interpolating
        Fth2(ind2,:)  = fthm(ind1,:)*Am(:,:,ind1) + fphim(ind1,:)*Bm(:,:,ind1);
        Fphi2(ind2,:) = -fthm(ind1,:)*Bm(:,:,ind1) + fphim(ind1,:)*Am(:,:,ind1);      
    end
end 
% take inverse fft in phi
Fth2 = ifft(Fth2,[],1);
Fphi2 = ifft(Fphi2,[],1);
end 