function [mPsin] = mPlmsin(L,theta)
% Normalized Legendre polynomial variant, l = 1,...,L
%
% m/sin(theta) \tilde P_l^m(cos(theta)) = 1/2 sqrt((l-1/2)/(l+1/2)) ...
%           (sqrt((l-m)(l-m-1)) \tilde P_{l-1}^{m+1}(cos(theta)) ...
%            + sqrt((l+m)(l+m-1)) \tilde P_{l-1}^{m-1}(cos(theta)) 
%
% L:        Maximum harmonic order L
% theta:	Spherical angles in radians
%
% mPsin:  	Derivative of normalized Legende polynomial evaluated at theta
%               Dimensions: size(theta) x N
%               N = L^2 + 2*L
%
% Dependencies: lm2ind

N = L^2 + 2*L;
theta = theta(:)';
Npoints = length(theta);
P = zeros(N,Npoints);   % P_l^m
mPsin = P;              % m P_l^m/sin(theta)

% load P_l^m for m >= 0
for l=1:L,  
    indlm0 = l^2 + l;  % lm2ind(l,0)
    indll = l^2 + 2*l; % lm2ind(l,l)
    ind = indlm0:indll;
    P(ind,:) = legendre(l,cos(theta),'norm');
end

% m P_l^m/sin(theta)
mPsin(lm2ind(1,1),:) = sqrt(3/4)*ones(1,Npoints);  % (l,m) = (1,1)
for l=2:L,
    c1 = 0.5*sqrt((l+0.5)/(l-0.5));
    for m=1:l,
       indlm = l^2 + l + m;             % lm2ind(l,m)
       indlm1mm1 = (l-1)^2 + l-1 + m-1; % lm2ind(l-1,m-1)
       mPsin(indlm,:) = c1*sqrt((l+m-1)*(l+m))*P(indlm1mm1,:); 
       if l-m > 1
          indlm1mp1 = (l-1)^2 + l-1 + m+1; % lm2ind(l-1,m+1)
          mPsin(indlm,:) = mPsin(indlm,:) + c1*sqrt((l-m)*(l-m-1))*P(indlm1mp1,:);
       end
    end
end

% negative m
for l=1:L,  
for m=-l:-1,
    indlm = l^2 + l + m;    % lm2ind(l,m)
    indlmm = l^2 + l - m;   % lm2ind(l,-m)
    mPsin(indlm,:) = -(-1)^(m)*mPsin(indlmm,:);
end
end

% extra (-1)^m
for l=1:L,
for m=-l:l,
    indlm = l^2 + l + m;    % lm2ind(l,m)
    mPsin(indlm,:) = (-1)^(m)*mPsin(indlm,:);
end
end
