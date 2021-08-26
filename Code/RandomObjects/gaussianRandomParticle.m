function [R] = gaussianRandomParticle(Theta,Phi,a,rmsh,Gamma_deg,L)
% Compute surface radii of a Gaussian random particle
%
% Theta,Phi:    [Any size] Spherical angles (radians) at which to compute the surface
% a:            Mean radius of surface (length)
% rmsh:         log-normal RMS of the surface (length)
% Gamma_deg:    Correlation angle (degrees)
% L:            [optional] Maximum harmonic degree L, defaults to heuristic
%
% R:            [Same size as Theta,Phi] Surface radius at (Thtea,Phi)
%
% Dependencies: coefclp, prodexpil
 
sigma = rmsh/a;                 % sigma parameter
beta = sqrt(log(sigma^2 + 1));  % beta parameter (std of zero mean Gaussian)
if nargin == 5
    L = round(275/Gamma_deg + 2.5); % number of harmonics given Gamma from heuristic
end

% harmonic indexing
tot = sum(0:L) + L + 1;         % total number of harmonics
lm = zeros(tot,2);              % store (l,m) pairs in lm
cnt = 0;
for l=0:L,
for m=0:l,
    cnt = cnt+1;
    lm(cnt,1) = l;
    lm(cnt,2) = m;
end
end

% pre-compute normalized associated Legendre polynomials
Leg = [];
for l=0:L,
    tmp = legendre(l,cos(Theta(:)),'norm');
    Leg = [Leg; tmp];
end

sz = size(Theta);
Theta = Theta(:)';
Phi = Phi(:)';
% make m index and phi gridded
[M P2] = ndgrid(0:L,Phi);
% pre-compute cosine and sine functions
C = cos(M.*P2);
S = sin(M.*P2);

% random draws for alm' and blm' from standard normal
% this indexing ensures that the seed counts up from the lowest harmonic 
tmp = randn(2*tot,1);
almp = tmp(1:2:end);
blmp = tmp(2:2:end);

% compute cl'
[clp] = coefclp(L,Gamma_deg);

% sum up the harmonics with correct coefficients
s = zeros(size(Theta));
ind = 0;
for l=0:L,
for m=0:l,
    ind = ind+1;
    if m==0
        const = sqrt(2);
    else
        const = 1;
    end
    tmp = clp(l+1)/const*(almp(ind)*C(m+1,:) + blmp(ind)*S(m+1,:));
    s = s + Leg(ind,:).*tmp;
end
end
s = beta*s;

% log-normal radius
R = a/sqrt(1 + sigma^2)*exp(s);
R = reshape(R,sz);

end

function [clp] = coefclp(L,Gamma)
% compute expansion coefficients of the correlation function
    lc = 2*sind(Gamma/2);
    x = 1/lc^2;
    tmp = prodexpil(L,x);
    clp = sqrt(2)*sqrt(2.*tmp);
end

function [S] = prodexpil(L,z)
% compute the product of exp(-z)*i_l(z) 
% by log transform with recursion
%
% S = exp(-z)*sqrt(pi/(2*z))*besseli((0:L)'+1/2,z)

lnz = log(z);
lnal0 = zeros(L+1,1);
% iterate a_{l+1,0)} inital conditions
lnal0(1) = -z;
for l = 0:(L-1),
   ind = l+1;
   lnal0(ind+1) = lnal0(ind) + lnz - log(2*l+3);    
end

% iterate the sums
S = zeros(L+1,1);
for l = 0:L,
    v = lnal0(l+1); 
    g = exp(v);
    first = v;
    k = 0;
    cnvg = 0;
    const = -log(2) + 2*lnz;
    while cnvg == 0 
        v = v + const - log(k+1) - log(2*(l+k)+3);
        g = g + exp(v);
        k = k + 1;
        if v < (first - 10)
            cnvg = 1;
        end
    end
    S(l+1) = g;
end
end