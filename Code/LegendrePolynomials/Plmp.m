function [plmp] = Plmp(L,x)
% Normalized associated Legendre polynomials derivative
% l = 0,...,L, m = -l:l.
%
%  d/dx \tildeP_l^m(x) = 1/(x^2-1)*(lx\tildeP_l^m(x) - ...
%               sqrt((l+1/2)/(l-1/2)*(l+m)(l-m))\tildeP_{l-1}^m(x)
%
% L:    maximum degree L
% x:    input arguments, x = [-1 1];
% str:  'scalar' to include (l,m) = (0,0)
%
% plmp:  Normalized associated Legendre polynomial derivative
%
% dependencies: Plm, lm2ind

x=x(:)';
plm = Plm(L,x);
plmp = zeros(size(plm));
for l = 0:L,
for m = -l:l,
    ind = lm2ind(l,m,'mono');
    if l == 0
        c2 = 0;
    elseif abs(m) <= l-1,        
        ind2 = lm2ind(l-1,m,'mono');
        c2 = -sqrt(l+1/2)/sqrt(l-1/2)*sqrt((l+m)*(l-m));
        c2 = c2.*plm(ind2,:);
    else
        c2 = 0;
    end
    c1 = l*x.*plm(ind,:);
    plmp(ind,:) = (1./(x.^2-1)).*(c1 + c2);
end
end


