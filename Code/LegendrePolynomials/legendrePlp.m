function [Plp] = legendrePlp(L,x)
% Legendre polynomial derivative computed with:
%
% P'_{l+1}(x) = (2l+1)P_l(x) + P'_{l-1}(x)
%
% L:    maximum degree L
% x:    evaluateion points -1 < x < 1
%
% Plp:  legendre polynomials evaluated for 0:L all x
%       Size (L+1) X length(x)
%
% Dependencies: legendrePl

if L < 0
    error('bad degree')
end
x = x(:)';
ind = or((x < -1),(x > 1));
if ~isempty(ind)
    x(ind) = NaN;
end
Plp = zeros(L+1,length(x));
if L == 0
    return
end
Plp(2,:) = 1;
if L == 1
    return
end
Pl = legendrePl(L,x);
for n=2:L,
   Plp(n+1,:) = (2*(n-1)+1)*Pl(n,:) + Plp(n-1,:);    
end


