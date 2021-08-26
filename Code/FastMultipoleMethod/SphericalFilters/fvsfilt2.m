function [Fth2 Fphi2] = fvsfilt2(Fth,Fphi,L,muj,wj,K,muk,wk)
% Fast vector spherical filter (interpolation or filtering)
% Based on fast scalar filter with correction terms.
% Assumes L <= K. Auto-detects filter or interplation based on size of Fth,
% Fphi. 
%
% L:        Lower degree harmonic (to interpolate up from, or filter down to) 
% K:        Upper degree harmonic (to interpolate up to, or filter down from)
%
% Interpolation
% Fth, Fphi:     input fields f(theta_j,phi_j):  (2L+1)x(L+1)
% Fth2, Fphi2:   output fields f(theta_k,phi_k): (2K+1)x(K+1)
%
% Filter
% Fth, Fphi:     input scalar field f(theta_k,phi_k):  (2K+1)x(K+1)
% Fth2, Fphi2:   output scalar field f(theta_j,phi_j): (2L+1)x(L+1)
%
% Dependences: fssfilt, vst, sphericalY, correctionTerms, h1

if L > K, error('error, L > K'), end
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
if M == I % interpolation 
    if N ~= J, error('bad array size'), end
    filt = 0;
    Fth2 = zeros(P,Q);
    Fphi2 = zeros(P,Q);
    mu1 = muj;
    w1 = wj;
    mu2 = muk;
    w2 = wk;
    I1 = I;
    I2 = P;
    J1 = J;
    J2 = Q;
elseif M == P; % filter 
    if N ~= Q, error('bad array size'), end
    filt = 1;
    Fth2 = zeros(I,J);
    Fphi2 = zeros(I,J);
	mu1 = muk;
    w1 = wk;
    mu2 = muj;
    w2 = wj;
    I1 = P;
    I2 = I;
    J1 = Q;
    J2 = J;
else
    error('bad array size')
end
sintheta1 = repmat(sqrt(1-(mu1(:)').^2),I1,1);
sintheta2 = repmat(sqrt(1-(mu2(:)').^2),I2,1);
phi2 = 2*pi*(0:(I2-1))/I2;
th2 = acos(mu2);
[Th Phi] = meshgrid(th2,phi2);
    
% fast scalar filter at L-1, fields are multiplied by sin of the current
% sampling
Fth_scalar = fssfilt(Fth.*sintheta1,L-1,muj,wj,K,muk,wk);
Fphi_scalar = fssfilt(Fphi.*sintheta1,L-1,muj,wj,K,muk,wk);

% compute expansion coefficients 
[blm, clm] = vst(Fth,Fphi,L,mu1,w1);

% indecies for L-1, L, L+1
Lp = L-1;
indLm1m = Lp^2 + Lp + (-Lp:Lp)';
Lp = L;
indLm = Lp^2 + Lp + (-Lp:Lp)';
Lp = L+1;
indLp1 = Lp^2 + Lp + (-Lp:Lp)';

% grab coefficents at L-1, L
bLm1m = blm(indLm1m);
bLm = blm(indLm);
cLm1m = clm(indLm1m);
cLm = clm(indLm);

% compute correction terms
[eL, eLp1, gL, gLp1] = correctionTerms(L,bLm1m,bLm,cLm1m,cLm);

% compute scalar spherical harmonics up to L+1 (will only use L and L+1)
ylm = sphericalY(L+1,Th,Phi);
grab = [indLm; indLp1];
Fth_corr = ylm(:,grab)*[eL; eLp1];
Fphi_corr = ylm(:,grab)*[gL; gLp1];
Fth_corr = reshape(Fth_corr,I2,J2);
Fphi_corr = reshape(Fphi_corr,I2,J2);

% sum the L-1 fast scalar filtered filed and the corrected field and divide
% by sin of the new sampling
Fth2 = (Fth_scalar + Fth_corr)./sintheta2;
Fphi2 = (Fphi_scalar + Fphi_corr)./sintheta2;
end

function [eL, eLp1, gL, gLp1] = correctionTerms(L,bLm1m,bLm,cLm1m,cLm)
% correction terms for K and K+1 coefficients
    tot1 = 2*L + 1;
    eL = zeros(tot1,1);
    gL = zeros(tot1,1);
    
    tot2 = 2*(L+1) + 1;
    eLp1 = zeros(tot2,1);
    gLp1 = zeros(tot2,1);
    
    % -(L-1):(L-1)
    m = ((-(L-1)):(L-1))';
    ind = (2:(tot1-1))';
    eL(ind) = bLm1m.*h1(L-1,m) + cLm(ind).*(1i*m);
	gL(ind) = -cLm1m.*h1(L-1,m) + bLm(ind).*(1i*m);
    eLp1(ind+1) = bLm(ind).*h1(L,m);
    gLp1(ind+1) = -cLm(ind).*h1(L,m); 
    
    % -L,L
    m = [-L, L]';
    ind = [1 tot1]';
    ind2 = [2 (tot2-1)]';
    eL(ind) = cLm(ind).*(1i*m);
    gL(ind) = bLm(ind).*(1i*m);
    eLp1(ind2) = bLm(ind).*h1(L,m);
    gLp1(ind2) = -cLm(ind).*h1(L,m);
end

% helper function
function h = h1(l,m)
    h = sqrt((l+0.5).*(l+1+m)./(l+1.5)./(l+1-m)).*l.*(l-m+1)./(2*l+1);
end
