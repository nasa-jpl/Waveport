function [ind] = lm2ind(l,m,str)
% Returns linear index for anglular harmonic number (l,m)
%
% l,m:  angular harmonic number
% str:  [optional] 'mono' for scalar wave function indexing
%
% ind:  linear index

if nargin == 3 
    if ~strcmp(str,'mono')
        error('bad string')
    end
    if sum(or(l<0,abs(m)>l))>0
        error('bad index')
    end
    ind = l.^2 + l + m + 1;
else
    if sum(or(l<1,abs(m)>l))>0
        error('bad index')
    end
    ind = l.^2 + l + m;
end


