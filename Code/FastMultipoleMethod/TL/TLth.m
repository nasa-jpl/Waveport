function [tl] = TLth(L,kX,costheta)
% FMM Far-field translation operator
%
% T_L = sum_{l=0}^L i^l(2l+1)h_l^{(1)}(kX)P_l(cos(theta))
%
% L:            Maximum order L
% kX:           Product of wavenumber k, and translation magnitude X
% costheta:     Dot product of khat and Xhat
%
% tl:           Translation operator values, size(costheta)
%
% Dependencies: sbesselh

sz = size(costheta);
costheta = costheta(:);
l = (0:L)';
H = sbesselh(l,kX);
H = ((1i).^(l).*(2*l+1)).*H;
tl = H(1)*ones(size(costheta));
if L == 0
    tl = reshape(tl,sz);
    return
end
tl = tl + H(2)*costheta;
if L == 1
    tl = reshape(tl,sz);
    return
end
Pl = 1.5*(costheta.^2) - 0.5; 
tl = tl + H(3)*Pl;
if L == 2
    tl = reshape(tl,sz);
    return
end
Plm1 = costheta;
for n=3:L,   
    Plm2 = Plm1;
    Plm1 = Pl;
    c1 = (2*n-1)/n;
    c2 = -(n-1)/n;
    Pl = c1*(costheta.*Plm1) + c2*Plm2; 
    tl = tl + H(n+1)*Pl;
end
tl = reshape(tl,sz);