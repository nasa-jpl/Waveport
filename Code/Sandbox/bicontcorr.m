function [G Calp Cs] = bicontcorr(r,lave,b,fv)
% Bicontinuous random media correlation and covariance functions
%
% r         Array of radial positions
% lave      Average length scale of heterogeneity
% b         Shape parameter
% fv        Volume fraction between 0 and 1
%
% G         Bicontinuous correlation function, \Gamma(r)
% Calp      Autocovariance function, C_{\alpha}(r)
% Cs        Autocvariance function, C_s(r)
%
% Dependencies: randsphere,sph2cart

% scale and shape parameters
sz = size(r);
r = r(:);
Nr = length(r);
kave = 2*pi/lave;
alp = erfinv(1-2*fv);

% compute C_s(r)
v = atan(kave*r/(b+1));
Cs = cos(v).^(b+1).*(v./sin(v)).*(sin(b*v)./(b*v));
Cs(r == 0) = 1;
lnCs = log(Cs);
signCs = sign(Cs);



% compute C_{\alpha}(r)
%Calp = zeros(size(Cs));
lnCm = -2*alp^2 - log(pi) - log(2);
M=100;
Calp = zeros(Nr,M);
Calp(:,1) = exp(lnCm)*Cs;
for m=2:M
    if m == 2
        %Hm = 2*alp;	% H1 = 2x
        %Hmm = 1;	% H0 = 1    
        hm = 2*alp; % H1/H0
    else
        %Hmp = 2*alp*Hm - 2*(m-2)*Hmm;
        %Hmm = Hm;
        %Hm = Hmp;
        hm = 2*alp - 2*(m-2)/hm;
    end
    lnCm = lnCm + 2*(log(abs(hm))) - log(m) - log(2);
%    Calp(:,m) = Calp(:,m-1) + exp(lnCm)*Cs.^m;
    Calp(:,m) = Calp(:,m-1) + signCs.^m.*exp(lnCm + m*lnCs);
%    Calp = Calp + signCs.^m.*exp(lnCm + m*lnCs);
end

% Calp2 = zeros(Nr,1);
% for n=1:Nr,
%     Calp2(n) = levint(Calp(n,:),M-2,1);
% end


Calp = Calp(:,end);

% % apply Shanks transform to the partial sums

% CalpST = zeros(size(Calp));
% grab1 = 1:(M-2);
% grab2 = 2:(M-1);
% grab3 = 3:M;
% D1 = Calp(:,grab3)-Calp(:,grab2);
% D2 = Calp(:,grab2)-Calp(:,grab1);
% CalpST = Calp(:,grab3) - (D1).^2./(D1 - D2);
% Calp = CalpST(:,end);

Calp(r==0) = fv-fv^2;

% correlation function
G = fv^2 + Calp;

