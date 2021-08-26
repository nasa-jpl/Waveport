function [plm] = Plm(L,x)
% Normalized associated Legendre polynomials evaluated at 
% l = 0,...,L, m = -l:l.
%
%  \tildeP_l^m(x) = (-1)^m sqrt((l+1/2)(l-m)!/(l+m)!) P_l^m(x)
%
% L:    maximum degree L
% x:    input arguments, x = [-1 1];
%
% plm:  Normalized associated Legendre polynomial
%
% dependencies: lm2ind

x=x(:);
N = length(x);
tot = L^2 + 2*L + 1;
plm = zeros(tot,N);
for l=0:L,
    indlm0 = l^2 + l + 1;   % lm2ind(l,0,'mono')
    indll = l^2 + 2*l + 1;  % lm2ind(l,l,'mono')
    ind = indlm0:indll;
    plm(ind,:) = legendre(l,x,'norm');    
end
for l=1:L,
for m=1:l,
   indlm = l^2 + l + m + 1;  % lm2ind(l,m,'mono')
   indlmm = l^2 + l - m + 1; % lm2ind(l,-m,'mono')
   plm(indlmm,:) = (-1)^m*plm(indlm,:);
end
end
for l=1:L,
for m=-l:l,
   indlm = l^2 + l + m + 1;  % lm2ind(l,m,'mono')
   plm(indlm,:) = (-1)^m*plm(indlm,:);
end
end
