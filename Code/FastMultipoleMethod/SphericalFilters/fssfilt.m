function [f2] = fssfilt(f,L,muj,wj,K,muk,wk)
% Fast scalar spherical filter, basic implementation 
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
% Dependencies: Plm, Plmp

if L > K
    disp('error, L > K')
    return
end
if L == K
    f2 = f;
    return
end
[M, N] = size(f);
J = length(muj);
I = 2*(J-1) + 1;
Q = length(muk);
P = 2*(Q-1) + 1;
if M == I % interpolation mode
    if N ~= J
        disp('bad array size')
        return
    end
    filt = 0;
    const = P/I;
    fmthp = zeros(P,Q);
    ind_sing = (Q+1)/2;
    sing_lim = J;
elseif M == P % filter mode
    if N ~= Q
        disp('bad array size')
        return
    end
    filt = 1;
    const = I/P;
    fmthp = zeros(I,J);
    % swap the quadrature nodes if filtering
    tmpmuj = muj;
    tmpwj = wj;
    muj = muk;
    wj = wk;
    muk = tmpmuj;
    wk = tmpwj;
    ind_sing = (J+1)/2;
    sing_lim = Q;
else
    disp('bad array size')
    return
end

% compute Legendre polynomials
PLj = Plm(L+1,muj);
PLk = Plm(L+1,muk);
PLjp = Plmp(L+1,muj);

% compute 1/(\mu_k - \mu_j) matrix
[Uk, Uj] = meshgrid(muk,muj);
M = 1./(Uk - Uj);

% if P and Q are odd, take care of singularity later
singularity = mod(J,2) && mod(Q,2);

% phi transform
fmth = fft(f,[],1);

% Combined forward/inverse Legendre transforms
indLmbase = L^2 + L + 1;
indLp1mbase = (L+1)^2 + (L+1) + 1;
for m=-L:L,
    indLm = indLmbase + m;
    indLp1m = indLp1mbase + m;
    if m >= 0
        indk = m + 1;
        indj = m + 1;
    else
        if ~filt
            indk = P + m + 1;
            indj = I + m + 1;
        else
            indk = I + m + 1;
            indj = P + m + 1;
        end
    end
        
    % terms of the sum
    bj = fmth(indj,:).*(wj).*PLj(indLm,:);    
    bjp1 = fmth(indj,:).*(wj).*PLj(indLp1m,:);
        
    % matrix-vector multiply
    coefj = bj*M;  
    coefjp1 = bjp1*M;
        
    % recompute one sum when the singularity exists       
    if singularity
        tmp1 = 0;
        tmp2 = 0;
        for n=1:sing_lim,
            c1 = fmth(indj,n)*wj(n);
            if n == (sing_lim+1)/2;
                tmp1 = tmp1 - c1*PLjp(indLm,n);
                tmp2 = tmp2 - c1*PLjp(indLp1m,n);
            else
                tmp1 = tmp1 + c1*PLj(indLm,n)/(muk(ind_sing)-muj(n));
                tmp2 = tmp2 + c1*PLj(indLp1m,n)/(muk(ind_sing)-muj(n));
            end
        end
        coefj(ind_sing) = tmp1;
        coefjp1(ind_sing) = tmp2;
    end
    bk = PLk(indLm,:);
    bkp1 = PLk(indLp1m,:);
    epsilon = sqrt(((L+1)^2 - m^2)/(4*(L+1)^2 - 1));
    fmthp(indk,:) = epsilon*(bkp1.*coefj - bk.*coefjp1);   
end

% inverse phi transform
f2 = const*ifft(fmthp,[],1);


