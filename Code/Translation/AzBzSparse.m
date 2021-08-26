function [row col A B] = AzBzSparse(L,Lp,k,rji,rgstr,normstr)
% Vector axial translation matrix
%
% L,Lp:     Maximum row, column harmonic 
% k:        Background wavenumber
% rji:      Translation distance along +z
% rgstr:    'rg' for beta = Rg{alpha}
% normstr:  'norm' for fully normalized wave functions
%
% row,col:  row and column indices of sparse matrix entries
% A, B:     sparse matrix entires
%
% Dependencies: lm2ind, sbesselh, sbesselj,
%               alphazSparse, NAzBz, indAzBz

if nargin == 4
    rgstr = [];
    normstr = [];
end
% sprase scalar axial translation matrix
[~, ~, alpz] = alphazSparse(L+1,Lp,k,rji,rgstr);
kr = k*rji;
A = zeros(NAzBz(L,Lp),1);
B = A;
row = A;
col = A; 

for n=1:Lp,
for l=1:L,
    lim = min(l,n);
    for m=-lim:lim,
        c1 = sqrt((l+m+1)*(l-m+1)/(2*l+1)/(2*l+3))/(l+1);
        c3 = m/(l*(l+1));
        if abs(m) <= l-1
            c2 = sqrt((l+m)*(l-m)/(2*l-1)/(2*l+1))/l;
            c2 = c2*alpz(indAlphaz(l-1,n,m,Lp));
        else
            c2 = 0;
        end
        inds = indAzBz(l,n,m,Lp);
        row(inds) = lm2ind(l,m);
        col(inds) = lm2ind(n,m);        
        in1 = indAlphaz(l+1,n,m,Lp);
        in2 = indAlphaz(l,n,m,Lp);        
        A(inds) = kr*(c1*alpz(in1) + c2) + alpz(in2);
        B(inds) = 1i*kr*c3*alpz(in2);        
    end
end
end
if nargin == 6 && strcmp(normstr,'norm')
    for n=1:Lp,
    for l=1:L,
        lim = min(l,n);
        for m=-lim:lim,
            const = sqrt(l*(l+1))/sqrt(n*(n+1));
            inds = indAzBz(l,n,m,Lp);
            A(inds) = const*A(inds);
            B(inds) = const*B(inds);
        end
    end
    end
end