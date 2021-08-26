function [tab] = lmtable(L,str)
% Table of indices for vector spherical harmonics
%
% L:    Maximum harmonic l
% str:  [optional] 'mono' for scalar wave function indexing
%
% tab:  2 column array with (l,m)
%       (l,m) = (column 1, column 2) 
     
if L < 0
    error('bad index')
end
tot = L^2 + 2*L;
startl = 1;
if nargin == 2
    if ~strcmp(str,'mono')
        error('bad string')
    end
    tot = tot + 1;
    startl = 0;
end
tab = zeros(tot,2);
cnt = 1;
for l=startl:L,
for m=(-l):l,
    tab(cnt,1) = l;
    tab(cnt,2) = m;
    cnt = cnt+1;
end
end