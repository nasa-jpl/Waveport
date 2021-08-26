function [A B] = AzBz(L,Lp,k,rji,rgstr,normstr)
% Vector axial translation matrix
%
% L,Lp:     Maximum row/column harmonic degree 
% k:        Background wavenumber
% rji:      Translation distance along +z
% rgstr:    'rg' for beta = Rg{alpha}
% normstr:  'norm' for fully normalized wave functions
%
% A,B:      Translation matrices size N1 x N2
%           N1 = L^2 + 2*L
%           N2 = Lp^2 + 2*Lp
%
% Dependencies: lm2ind, alphaz

if nargin == 4
    rgstr = [];
    normstr = [];
end
% compute scalar axial translation
alpz = alphaz(L+1,Lp,k,rji,rgstr);

tot1 = L^2 + 2*L;
tot2 = Lp^2 + 2*Lp;
kr = k*rji;
A = zeros(tot1,tot2);
B = A;
for n=1:Lp,
for l=1:L,
    lim = min(l,n);
    for m=-lim:lim,
        lm = lm2ind(l,m,'mono');
        lp1m = lm2ind(l+1,m,'mono');
        nm = lm2ind(n,m,'mono');
        c1 = sqrt((l+m+1)*(l-m+1)/(2*l+1)/(2*l+3))/(l+1);
        c3 = m/(l*(l+1));
        if abs(m) <= l-1
            lmm1 = lm2ind(l-1,m,'mono');
            c2 = sqrt((l+m)*(l-m)/(2*l-1)/(2*l+1))/l;
            c2 = c2*alpz(lmm1,nm);
        else
            c2 = 0;
        end
        lm2 = lm2ind(l,m);
        nm2 = lm2ind(n,m);
        A(lm2,nm2) = kr*(c1*alpz(lp1m,nm) + c2) + alpz(lm,nm);
        B(lm2,nm2) = 1i*kr*c3*alpz(lm,nm);
    end
end
end
if nargin == 6 && strcmp(normstr,'norm')
    for n=1:Lp,
    for l=1:L,
        lim = min(l,n);
        for m=-lim:lim,
            const = sqrt(l*(l+1))/sqrt(n*(n+1));
            A(lm2ind(l,m),lm2ind(n,m)) = const*A(lm2ind(l,m),lm2ind(n,m));
            B(lm2ind(l,m),lm2ind(n,m)) = const*B(lm2ind(l,m),lm2ind(n,m));
        end
    end
    end
end