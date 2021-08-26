function [l, m] = ind2lm(ind,str)
% Returns angular harmonic number (l,m) from linear index
%
% ind:  linear index
% str:  [optional] 'mono' for scalar wave function indexing
%
% l,m:  angular harmonic number

if sum(ind<1)>0
    error('bad index')
end
if nargin == 2 
    if ~strcmp(str,'mono')
        error('bad string')
    end
    ind = ind - 1;
end
l = floor(sqrt(ind));
m = ind - l.^2 - l;


