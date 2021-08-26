function [row col alpha] = alphazSparse(L,Lp,k,rji,rgstr)
% Scalar axial translation matrix computed as sparse matrix
%
% L,Lp:     Maximum row, column degree harmonic 
% k:        Background wavenumber
% rji:      Translation distance along +z
% rgstr:    'rg' for beta = Rg{alpha}, uses spherical Bessel j_l(kr)
%
% row,col:  row and column indices of sparse matrix entries
% alpha:    sparse matrix entries
%
% Dependencies: lm2ind, sbesselh, sbesselj,
%                   Nalphaz, Nalphazfill, indAlphazFill, indAlphaz

% Determine maximum degree and transpose condition
flip = 0;
if Lp > L
    tmp = Lp;
    Lp = L;
    L = tmp;
    flip = 1;
end
L3 = L + Lp;
str = 'mono';
NLLP = Nalphaz(L,Lp);
NN = NLLP + Nalphazfill(L,Lp);
alpha = zeros(NN,1);
row = zeros(NN,1);
col = zeros(NN,1);

% Initialize, alpha_{l0;00}
nk = lm2ind(0,0,str);
for l=0:L3,
   lm = lm2ind(l,0,str);
    if nargin == 5 && strcmp(rgstr,'rg')
        bess = sbesselj(l,k*rji);
    else
        bess = sbesselh(l,k*rji);
    end
    inds = indAlphazFill(l,0,0,L,Lp);
    row(inds) = lm;
    col(inds) = nk;
    alpha(inds) = (-1)^l*sqrt(2*l+1)*bess;
end

% alpha_{l,n+1;n+1,n+1}
for n = 0:(Lp-1),
for l = (n+1):(L3-n-1),
    lnp1 = lm2ind(l,n+1,str);
    np1np1 = lm2ind(n+1,n+1,str);
    if l-1 >= 0
        c2 = sqrt((l+n+1)*(l+n)/(2*l-1)/(2*l+1));
        c2 = c2*alpha(indAlphazFill(l-1,n,n,L,Lp));
    else
        c2 = 0;
    end
    c1 = sqrt((2*n+3)/(2*(n+1)));
    c3 = sqrt((l-n+1)*(l-n)/(2*l+3)/(2*l+1));
    inds = indAlphazFill(l,n+1,n+1,L,Lp);
    row(inds) = lnp1;
    col(inds) = np1np1;
    alpha(inds) = c1*(c2 + c3*alpha(indAlphazFill(l+1,n,n,L,Lp)));
end    
end

% alpha_{l,-n;n,-n}
for n = 1:Lp,
for l = n:(L3-n),
    lmn = lm2ind(l,-n,str);
    nmn = lm2ind(n,-n,str);
    inds = indAlphazFill(l,n,-n,L,Lp);
    row(inds) = lmn;
    col(inds) = nmn;
    alpha(inds) = alpha(indAlphazFill(l,n,n,L,Lp));
end
end

% alpha_{lm,n+1,m}, m = +/-n
for n=0:(Lp-1),
for l=0:(L3-n-1),
    lim = min(l,n);
    for m=[-lim,lim],
        lm = lm2ind(l,m,str);
        np1m = lm2ind(n+1,m,str);
        c1 = sqrt(2*n+3);
        if abs(m) <= l-1
            c2 = sqrt((l+n)*(l-n)/(2*l-1)/(2*l+1));
            c2 = c2*alpha(indAlphazFill(l-1,n,m,L,Lp));
        else
            c2 = 0;
        end
        c3 = sqrt((l+n+1)*(l-n+1)/(2*l+3)/(2*l+1));
        inds = indAlphazFill(l,n+1,m,L,Lp);
        row(inds) = lm;
        col(inds) = np1m;
        alpha(inds) = c1*(c2 - c3*alpha(indAlphazFill(l+1,n,m,L,Lp)));
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
        c0 = sqrt((n+m+1)*(n-m+1)/(2*n+1)/(2*n+3));
        c1 = sqrt((n+m)*(n-m)/(2*n+1)/(2*n-1));
        if abs(m) <= l-1
            c2 = sqrt((l+m)*(l-m)/(2*l+1)/(2*l-1));
            c2 = c2*alpha(indAlphazFill(l-1,n,m,L,Lp));
        else
            c2 = 0;
        end
        c3 = sqrt((l+m+1)*(l-m+1)/(2*l+3)/(2*l+1));
        inds = indAlphazFill(l,n+1,m,L,Lp);
        row(inds) = lm;
        col(inds) = np1m;        
        alpha(inds) = (1/c0)*(c1*alpha(indAlphazFill(l,n-1,m,L,Lp)) + c2 ...
                        - c3*alpha(indAlphazFill(l+1,n,m,L,Lp)));    
    end   
end
end

% trim the array
alpha = alpha(1:NLLP);
row = row(1:NLLP);
col = col(1:NLLP);

% if transpose, apply (-1)^(l+n) and reorder
% to be consistent with row-major indexing from indAlphaz
if flip
    rtmp = zeros(NLLP,1);
    ctmp = rtmp;
    atmp = rtmp;
    for n=0:Lp,
    for l=0:L,
        lim = min(l,n);
        for m=-lim:lim,
            ind1 = indAlphaz(l,n,m,Lp);
            ind2 = indAlphaz(n,l,m,L);
            atmp(ind2) = (-1)^(n+l)*alpha(ind1);
            ctmp(ind2) = row(ind1);
            rtmp(ind2) = col(ind1);
        end
    end
    end
    alpha = atmp;
    row = rtmp;
    col = ctmp;    
end

   

