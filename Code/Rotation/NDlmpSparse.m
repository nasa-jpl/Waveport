function [Nl] = NDlmpSparse(L,str)
% Number of non-zero matrix entries for Dlmp rotation matrix
%
% L:    Maximum harmonic L
% str:  'mono' for monopole
%
% Nl:   Number of non-zero matrix entries

if nargin == 2 && strcmp(str,'mono')
    if L < 0, error('bad index'), end    
    Nl = round(4/3*L^3 + 4*L^2 + 11/3*L) + 1;
else
    if L < 1, error('bad index'), end
    Nl = round(4/3*L^3 + 4*L^2 + 11/3*L);
end

