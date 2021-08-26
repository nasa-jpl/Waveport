function [Pl] = legendrePl(L,x)
% Legendre polynomial computed with:
%
% (l+1)P_{l+1}(x) = (2l+1)xP_l(x) - lP_{l-1}(x)
%
% L:    maximum degree L
% x:    evaluateion points -1 < x < 1
%
% Pl:   legendre polynomials evaluated for 0:L all x
%       Size (L+1) X length(x)

if L < 0
    error('bad degree')
end
x = x(:)';
ind = or((x < -1),(x > 1));
if ~isempty(ind)
    x(ind) = NaN;
end
Pl = zeros(L+1,length(x));
Pl(1,:) = 1;
if L == 0
    return
end
Pl(2,:) = x;
if L == 1
    return
end
for n=1:(L-1),
   c1 = (2*n+1)/(n+1);
   c2 = -n/(n+1);
   Pl(n+2,:) = c1*(x.*Pl(n+1,:)) + c2*Pl(n,:);    
end

