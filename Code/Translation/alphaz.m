function [alpha] = alphaz(L,Lp,k,rji,rgstr)
% Scalar axial translation matrix
%
% L,Lp:     Maximum row, column degree harmonic 
% k:        Background wavenumber
% rji:      Translation distance along +z
% rgstr:    'rg' for beta = Rg{alpha}, uses spherical Bessel j_l(kr)
%
% alpha:    Translation matrix of size N1xN2
%           N1 = L^2 + 2*L + 1
%           N2 = Lp^2 + 2*Lp + 1
%
% Dependencies: lm2ind, sbesselh, sbesselj

% Determine maximum order and transpose condition
flip = 0;
if Lp > L
    tmp = Lp;
    Lp = L;
    L = tmp;
    flip = 1;
end
L3 = L + Lp;
tot1 = L^2 + 2*L + 1;
tot2 = Lp^2 + 2*Lp + 1;
tot3 = L3^2 + 2*L3 + 1;
alpha = zeros(tot3,tot2);
str = 'mono';

% Initialize, alpha_{l0;00}
nk = lm2ind(0,0,str);
for l=0:L3,
    lm = lm2ind(l,0,str);
    if nargin == 5 && strcmp(rgstr,'rg')
        bess = sbesselj(l,k*rji);
    else
        bess = sbesselh(l,k*rji);
    end
    alpha(lm,nk) = (-1)^l*sqrt(2*l+1)*bess;
end

% alpha_{l,n+1;n+1,n+1}
for n = 0:(Lp-1),
for l = (n+1):(L3-n-1),
    lnp1 = lm2ind(l,n+1,str);
    np1np1 = lm2ind(n+1,n+1,str);
    nn = lm2ind(n,n,str);
    lp1n = lm2ind(l+1,n,str);
    if l-1 >= 0
        lm1n = lm2ind(l-1,n,str);
        c2 = sqrt((l+n+1)*(l+n)/(2*l-1)/(2*l+1));
        c2 = c2*alpha(lm1n,nn);
    else
        c2 = 0;
    end
    c1 = sqrt((2*n+3)/(2*(n+1)));
    c3 = sqrt((l-n+1)*(l-n)/(2*l+3)/(2*l+1));
    alpha(lnp1,np1np1) = c1*(c2 + c3*alpha(lp1n,nn));
end    
end

% alpha_{l,-n;-n,-n}
for n = 1:Lp,
for l = n:(L3-n),
    ln = lm2ind(l,n,str);
    nn = lm2ind(n,n,str);
    lmn = lm2ind(l,-n,str);
    nmn = lm2ind(n,-n,str);
    alpha(lmn,nmn) = alpha(ln,nn); 
end
end

% alpha_{lm,n+1,m}, m = +/-n
for n=0:(Lp-1),
for l=0:(L3-n-1),
    lim = min(l,n);
    for m=[-lim,lim],
        lm = lm2ind(l,m,str);
        np1m = lm2ind(n+1,m,str);
        nm = lm2ind(n,m,str);
        lp1m = lm2ind(l+1,m,str);
        c1 = sqrt(2*n+3);
        if abs(m) <= l-1
            lm1m = lm2ind(l-1,m,str);
            c2 = sqrt((l+n)*(l-n)/(2*l-1)/(2*l+1));
            c2 = c2*alpha(lm1m,nm);
        else
            c2 = 0;
        end
        c3 = sqrt((l+n+1)*(l-n+1)/(2*l+3)/(2*l+1));
        alpha(lm,np1m) = c1*(c2 - c3*alpha(lp1m,nm));
    end   
end
end

% alpha_{lm,n+1,m}, m ~= +/-n
for n=1:(Lp-1),
for l=0:(L3-n-1), 
    lim = min(l,n-1);
    for m=-lim:lim,
        lm = lm2ind(l,m,str);
        np1m = lm2ind(n+1,m,str);
        nm1m = lm2ind(n-1,m,str);
        nm = lm2ind(n,m,str);
        lp1m = lm2ind(l+1,m,str);
        c0 = sqrt((n+m+1)*(n-m+1)/(2*n+1)/(2*n+3));
        c1 = sqrt((n+m)*(n-m)/(2*n+1)/(2*n-1));
        if abs(m) <= l-1
            lm1m = lm2ind(l-1,m,str);
            c2 = sqrt((l+m)*(l-m)/(2*l+1)/(2*l-1));
            c2 = c2*alpha(lm1m,nm);
        else
            c2 = 0;
        end
        c3 = sqrt((l+m+1)*(l-m+1)/(2*l+3)/(2*l+1));
        alpha(lm,np1m) = (1/c0)*(c1*alpha(lm,nm1m) + c2 - c3*alpha(lp1m,nm));    
    end   
end
end

% trim the matrix
alpha = alpha(1:tot1,:);

% apply (-1)^(l+n) if needed
if flip
    for n=0:Lp,
    for l=0:L,
        if mod(n+l,2) == 1
            lim = min(l,n);
            for m=-lim:lim,
                lm = lm2ind(l,m,str);
                nm = lm2ind(n,m,str);
                alpha(lm,nm) = (-1)^(n+l)*alpha(lm,nm);
            end
        end
    end
    end
    alpha = alpha.';
end

