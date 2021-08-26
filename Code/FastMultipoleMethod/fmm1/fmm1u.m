function [u] = fmm1u(t,ti)
% Basis functions for 1D fast-multipole method
%
% u_j(t) = prod_{k=1, k~=j}^p (t-t_k)/(t_j-t_k)
%
% t:    input arguments
% ti:   p Chebychev coefficients on x = [-1 1]
%
% u:    basis function evaluated at t, size: length(t) x p

t = t(:);
p = length(ti);
n = length(t);
u = zeros(n,p);
for j=1:p
    tmp = ones(n,1);
    for k=1:p,
        if k~=j
            tmp = tmp.*((t-ti(k))/(ti(j)-ti(k)));
        end    
    end
    u(:,j) = tmp;
end


