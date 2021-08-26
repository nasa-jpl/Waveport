function [S] = fmm1prep(xk,xj,epsilon,s)
% Preparatory function for 1D fast multipole method routine fmm1
%
% xk:       source points, N x 1
% xj:       observation points, N x 1
% epsilon:  precision of the computation
% s:        [optional] number of points in each bin at finest level
%           default: s = 2*p
%
% Outputs stored in structure S (also ak and xj):
%
% ML,MR,SL,SR,T1,T2,T3,T4:
%           p x p translation matrices
% ti:       p x 1 Chebyshev coefficents on x = [-1 1]
% p:        maximum degree of expansions
% nlevs:    number of levels
% r:        radius of subintervals at finest level (half the width)
% cen:      box centers at finest level
% xj_ind:   box index of observation points at finest level
% xk_ind:   box index of source points at finest level
% Phi, Psi: p x B array, preallocated storage for expansion coefficients  
% B:        total number of boxes

% dependencies: fmm1u

p = ceil(-log(epsilon)/log(5));
if nargin == 3
    s = 2*p;
end
xj = xj(:);
xk = xk(:);
Nj = length(xj);
Nk = length(xk);
nlevs = ceil(log2(Nk/s));
if nlevs < 1 
   nlevs = 1;
end
B = 2^(nlevs+1) - 2;
Phi = zeros(p,B);
Psi = zeros(p,B);

% parameters at finest level
nbox = 2^(nlevs);
width = (1/2)^(nlevs-1);
r = width/2;
cen = (0:(nbox-1))*width - 1 + r;

% Chebyshev coefficients
ti = cos((2*(1:p)-1)*pi/p/2);

% box index at the finest level of the points
xj_ind = ceil((xj + 1)/2*nbox);
xj_ind(xj_ind == 0) = 1;
xk_ind = ceil((xk + 1)/2*nbox);
xk_ind(xk_ind == 0) = 1;

xj_ind_list = cell(nbox,1);
xk_ind_list = cell(nbox,1);
xj_notempty = [];
xk_notempty = [];
nn_ind_list = cell(nbox,1);
uj_local = cell(nbox,1);
uk_far = cell(nbox,1);
for i=1:nbox,
    xj_ind_list{i} = find(xj_ind == i); 
    xk_ind_list{i} = find(xk_ind == i);
    if ~isempty(xj_ind_list{i})
        xj_notempty = [xj_notempty; i];
        uj_local{i} = fmm1u((xj(xj_ind_list{i})-cen(i))/r,ti);
    end
    if ~isempty(xk_ind_list{i})
        xk_notempty = [xk_notempty; i];
        tmp = zeros(p,length(xk_ind_list{i}));    
        for pp=1:p,
            tmp(pp,:) = ti(pp)./(3*r-ti(pp)*(xk(xk_ind_list{i})-cen(i))); 
        end
        uk_far{i} = tmp;
    end    
    ind1 = find(xk_ind == i-1);
    ind2 = find(xk_ind == i);
    ind3 = find(xk_ind == i+1);
    indk = [ind1; ind2; ind3];
    nn_ind_list{i} = indk;
end

% compute contributions from neighbor boxes directly
% store as sparse matrix
row = [];
col = [];
elm = [];
for i=1:nbox,
    indj = xj_ind_list{i};
    indk = nn_ind_list{i};
    [indK, indJ] = ndgrid(indk,indj);
    M = 1./(xj(indJ) - xk(indK));
    M(isnan(M)) = 0;
    row = [row; indJ(:)];
    col = [col; indK(:)];
    elm = [elm; M(:)];
end
M = sparse(row,col,elm,Nj,Nk);

% M, S, T
ML = fmm1u(3*ti./(6+ti),ti);
MR = fmm1u(3*ti./(6-ti),ti);
SL = fmm1u((ti-1)/2,ti);
SR = fmm1u((ti+1)/2,ti);
T1 = fmm1u(3./(ti+6),ti);
T2 = fmm1u(3./(ti+4),ti);
T3 = fmm1u(3./(ti-4),ti);
T4 = fmm1u(3./(ti-6),ti); 

% store in structure
S.xk = xk;  S.xj = xj;  S.ML = ML;  S.MR = MR; S.SL = SL;
S.SR = SR;  S.T1 = T1;  S.T2 = T2;  S.T3 = T3;
S.T4 = T4;  S.ti = ti;  S.Nk = Nk;  S.Nj = Nj;
S.p = p; S.nlevs = nlevs; S.nbox = nbox; S.r = r;
S.cen = cen; S.xj_ind = xj_ind_list; S.xk_ind = xk_ind_list;
S.nn_list = nn_ind_list; S.uj_local = uj_local;
S.uk_far = uk_far; S.Phi = Phi; S.Psi = Psi; S.B = B; S.M = M; 
S.xj_notempty = xj_notempty; S.xk_notempty = xk_notempty;

end


