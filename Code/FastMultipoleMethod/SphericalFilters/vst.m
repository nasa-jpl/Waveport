function [blm, clm] = vst(Fth,Fphi,L,muj,wj,normstr,plm,dplm)
% Forward vector spectral transform
% Defaults to partially normalized vector spherical hamornics
%
% L:        maximum degree L<=L'
% Fth:      IxJ sampled theta component F_theta(theta,phi)
% Fphi:     IxJ sampled phi component F_phi(theta,phi)
%           I = 2L'+1, J = L'+1
%           phi_i = 2*pi*(0:(I-1))/I
%           theta = arccos(mu_j)
% muj:      J Gaussian quadrature nodes
% wj:       J Gaussian quadrature weights
% normstr:  [optional] 'norm' for fully normalized spherical harmonics
%                      'none' for no normalization
% plm,dplm: [optional] precomputed Legendre polynomials
%           plm = Plm(L,muj);
%           dplm = Plmp(L,muj); 
%
% blm,clm:  L^2+2*L spectral coefficients, linearly indexed
%
% dependencies: Plm, Plmp

[I, J] = size(Fth);
tot = L^2 + 2*L;

% fft in \phi
fthm = sqrt(2*pi)/I*fft(Fth,[],1);
fphim = sqrt(2*pi)/I*fft(Fphi,[],1);

% Legendre polynomials and derivative
if nargin <= 6
    plm = Plm(L,muj); % includes monopole
    dplm = Plmp(L,muj); % includes monopole
elseif nargin == 8 && ~exist('plm') && ~exist('dplm')
    disp('bad inputs for precomputation')
    return
end

muj = muj(:)';
wj = wj(:)';
sq = sqrt(1-muj.^2);
blm = zeros(tot,1);
clm = zeros(tot,1);
for l=1:L,
for m=-l:l,
    ind = l^2 + l + m; 
    c1 = -sq.*dplm(ind+1,:); % monopole indexing  
    c2 = -1i*m*plm(ind+1,:)./sq; % monopole indexing  
    if m >= 0
        ind2 = m + 1;
    else
        ind2 = I + m + 1;
    end
    if nargin >= 6 && strcmp(normstr,'none')
        const = 1;
    elseif nargin >= 6 && strcmp(normstr,'norm')
        const = 1/sqrt(l*(l+1));
    else
        const = 1/(l*(l+1));
    end
    blm(ind) = const*sum((c1.*fthm(ind2,:) + c2.*fphim(ind2,:)).*wj);
    clm(ind) = const*sum((c2.*fthm(ind2,:) - c1.*fphim(ind2,:)).*wj);
end 
end


