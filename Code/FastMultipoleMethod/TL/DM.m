function [dm] = DM(M,theta)
% Dirichlet kernel for FMM translation operator interpolation
%
% D_M(theta) = sin((2M+1)theta/2)((2M+1)sin(theta/2))
%
% M:        number of sample points between (0,pi), not including ends
% theta:    interpolation point, angle in radians
%
% dm:       Dirichlet kernel 

tmp = 2*M + 1;
thover2 = theta/2;
dm = sin(tmp*thover2)./(tmp*sin(thover2));
ind = find(theta == 0);
if ~isempty(ind)
    dm(ind) = 1;
end