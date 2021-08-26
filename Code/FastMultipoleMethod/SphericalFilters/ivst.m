function [Fth,Fphi] = ivst(blm,clm,L,muj,normstr,plm,dplm)
% Inverse vector spectral transform
%
% L:        Maximum degree L<=L'
% blm,clm:  L^2 + 2*L spectral coefficients on B_{lm} and C_{lm}
% muj:      J Gaussian quadrature nodes, muj = cos(theta_j) (L' = J-1 harmonics)
% normstr:  [optional] 'norm' for fully normalized spherical harmonics
% plm,dplm: [optional] precomputed Legendre polynomials
%           plm = Plm(L,muj);
%           dplm = Plmp(L,muj);
%
% Fth:      I'xJ' sampled theta component F_theta(theta,phi)
% Fphi:     I'xJ' sampled phi component F_phi(theta,phi)
%           I' = 2L'+1, J = L'+1
%
% dependencies: Plm, Plmp,

J = length(muj);
I = 2*(J-1) + 1;
if nargin <= 5
    plm = Plm(L,muj); % includes monopole
    dplm = Plmp(L,muj); % includes monopole
elseif nargin == 7 && ~exist('plm') && ~exist('dplm')
    disp('bad inputs for precomputation')
    return
end

Fth = zeros(I,J);
Fphi = zeros(I,J);

sq = sqrt(1-muj.^2)';
for m = -L:L,
    if m >= 0
        ind2 = m + 1;
    else
        ind2 = I + m + 1;
    end    
    if m == 0
        lowerlim = 1;
    else
        lowerlim = abs(m);
    end
    for l = lowerlim:L,
        if nargin >= 5 && strcmp(normstr,'norm')
            const = 1/sqrt(l*(l+1));
        else
            const = 1;
        end
        ind = l^2 + l + m;
        c1 = -sq.*dplm(ind+1,:);
        c2 = 1i*m*plm(ind+1,:)./sq;
        Fth(ind2,:)   = Fth(ind2,:) + const*(blm(ind)*c1 + clm(ind)*c2);
        Fphi(ind2,:)  = Fphi(ind2,:) + const*(blm(ind)*c2 - clm(ind)*c1);     
    end            
end 

Fth = I/sqrt(2*pi)*ifft(Fth,[],1);
Fphi = I/sqrt(2*pi)*ifft(Fphi,[],1);


