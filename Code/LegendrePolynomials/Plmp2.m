function [dP] = Plmp2(L,theta)
% Derivative of normalized Legendre polynomial,
% l = 1,...,L
%
% d/dtheta \tilde P_l^m(cos(theta)) = ...
%           1/2(sqrt((l+m)(l-m+1)) \tilde P_l^{m-1}(cos(theta)) ...
%            + sqrt((l-m)(l+m+1)) \tilde P_l^{m+1}(cos(theta)) 
%
% L:      Maximum harmonic order L
% theta:	Spherical angles in radians
% 
% dP:   	Derivative of normalized Legende polynomial evaluated at theta
%         Dimensions: size(theta) x N
%             N = L^2 + 2*L
%
% Dependencies: lm2ind

N = L^2 + 2*L;
theta = theta(:)';
Npoints = length(theta);
P = zeros(N,Npoints);   % P_l^m
dP = P;                 % d/dtheta P_l^m

% load P_l^m for m >= 0
for l=1:L,  
    indlm0 = l^2 + l;  % lm2ind(l,0)
    indll = l^2 + 2*l; % lm2ind(l,l)
    ind = indlm0:indll;
    P(ind,:) = legendre(l,cos(theta),'norm');
end

% d/dtheta P_l^m
for l=1:L,
    indlm0 = l^2 + l;       % lm2ind(l,0)
    indlm1 = l^2 + l + 1;   % lm2ind(l,1)
    dP(indlm0,:) = -sqrt(l*(l+1))*P(indlm1,:); 
    for m=1:l,
        indlm = l^2 + l + m;  % lm2ind(l,m)
        if m-1 >= 0
            indlmm1 = l^2 + l + m - 1;
            dP(indlm,:) = 0.5*sqrt((l+m)*(l-m+1))*P(indlmm1,:);
        end
        if m+1 <= l
            indlmp1 = l^2 + l + m + 1; % lm2ind(l,m+1)
            dP(indlm,:) = dP(indlm,:) - 0.5*sqrt((l-m)*(l+m+1))*P(indlmp1,:);
        end
    end
end

% negative m
for l=1:L,
for m=-l:-1,
    indlm = l^2 + l + m; % lm2ind(l,m)
    indlmm = l^2 + l -m; % lm2ind(l,-m)
    dP(indlm,:) = (-1)^(m)*dP(indlmm,:);
end
end

% extra (-1)^m
for l=1:L,
for m=-l:l,
    indlm = l^2 + l + m; % lm2ind(l,m)
    dP(indlm,:) = (-1)^(m)*dP(indlm,:);
end
end

