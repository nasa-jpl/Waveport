function [cj pj bj] = ebcmbesselfunc(tot,L,l,l2,k,r,rgstr)

% kr
kr = k*r;

% storage
cj = zeros(1,tot);
pj = zeros(1,tot);
bj = zeros(1,tot);
cjtmp = zeros(1,L);
pjtmp = zeros(1,L);
bjtmp = zeros(1,L);

% evalute bessel functions
if strcmp(rgstr,'rg')
    cjtmp(l2) = sbesselj(l2,kr);      % j_l(kr), handles lim r->0 
    pjtmp(l2) = sbesselj2(l2,kr);     % j_l(kr)/kr, "
    bjtmp(l2) = sbesseljp(l2,kr);     % j'_l(kr), " 
else
    cjtmp(l2) = sbesselh(l2,kr);       % h_l^(1)(kr)
    pjtmp(l2) = sbesselh(l2,kr)./kr;	% h_l^(1)(kr)/kr
    bjtmp(l2) = sbesselhp(l2,kr);      % h'_l^(1)(kr)
end
% complete the B function
bjtmp = pjtmp + bjtmp;

% apply scale factor to P function
pjtmp = sqrt(l2.*(l2+1)).*pjtmp;

% map into the full array
cj(:) = cjtmp(l);
pj(:) = pjtmp(l);
bj(:) = bjtmp(l);
end


