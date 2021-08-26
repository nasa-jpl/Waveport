function [Tmm Tmn Tnm Tnn] = ebcm(L,X,Y,Z,dS,nx,ny,nz,k1,k2,mu1,mu2)
% T-matrix computed with extended boundary condition method
%
% L:        Maximum degree harmonic
% X,Y,Z:    [any size, same size] Cartesian points of surface discritization
% dS:       [any size, same as X,Y,Z] Differential surface area for each surface point
% nx,ny,nz: [any size, same as X,Y,Z] Outward surface normal unit vector
% k1,k2:    Wavenumbers in outer and inner regions
% mu1,mu2:  [optional] Permeability of outer and inner regions, realative or absolute, default 1
% 
% Tmm, Tmn, Tnm, Tnn:   [NxN] Block T-matrices, N = L^2 + 2L
%
% Dependencies: lmtable,cart2sph,sph2cart,BC,sphericalY,ebcmbessel,ebcmprod 

if nargin == 10
    mu1 = 1;
    mu2 = 1;
end

% total number of surface points and harmonics
N = numel(X);
tot = L^2 + 2*L;

% storage for Q matrices
Q = zeros(2*tot,2*tot);
RgQ = zeros(2*tot,2*tot);
ZZ = zeros(tot,tot);
Qmm = ZZ;
Qmn = ZZ;
Qnm = ZZ;
Qnn = ZZ;
RgQmm = ZZ;
RgQmn = ZZ;
RgQnm = ZZ;
RgQnn = ZZ;

% subindices of M and N blocks in Q and RgQ
indM = 1:tot;
indN = (tot+1):(2*tot);

% precompute harmonic linear index
[tab] = lmtable(L);
l = tab(:,1);
l2 = (1:L)';

% multiplying constant
const = mu1*k2/(mu2*k1);

% create n_hat for dot products later
n_hat = [nx(:) ny(:) nz(:)];

% loop over surface points (slow but saves memory)
for n=1:N,
    % convert surface point to spherical coordinates
    [r th ph] = cart2sph(X(n),Y(n),Z(n));  
    
    % convert spherical unit vectors at this point to Cartesian components
    [r_hat_x r_hat_y r_hat_z]      = sph2cart(r,th,ph,1,0,0);
    [th_hat_x th_hat_y th_hat_z]   = sph2cart(r,th,ph,0,1,0);
    [ph_hat_x ph_hat_y ph_hat_z]   = sph2cart(r,th,ph,0,0,1);

    % dot product of surface normal and spherical unit vectors at the same point
    % multiply by the differential surface element here
    n_dot_r  = dS(n)*dot(n_hat(n,:),[r_hat_x r_hat_y r_hat_z],2);
    n_dot_th = dS(n)*dot(n_hat(n,:),[th_hat_x th_hat_y th_hat_z],2);
    n_dot_ph = dS(n)*dot(n_hat(n,:),[ph_hat_x ph_hat_y ph_hat_z],2);

    % unconjugated vector spherical harmonics up to degree L and +/- m (row arrays)
    [V2, U2, ~, ~] = BC(L,th,ph,'norm');
    W2 = sphericalY(L,th,ph);

    % conjugated vector spherical harmonics (row arrays)
    V1 = conj(V2);
    U1 = conj(U2);
    W1 = conj(W2);
    
    % bessel functions (row arrays)
    [c1 p1 b1]          = ebcmbessel(tot,L,l,l2,k1,r,[]);
    [rgc1 rgp1 rgb1]    = ebcmbessel(tot,L,l,l2,k1,r,'rg');
    [rgc2 rgp2 rgb2]    = ebcmbessel(tot,L,l,l2,k2,r,'rg');
    
    % cross product pre-multiplications (column arrays for 1, row arrays for 2)
    [cU1 cV1 bU1 bV1 pW1]           = ebcmprod(c1.',p1.',b1.',U1.',V1.',W1.');
    [rgcU1 rgcV1 rgbU1 rgbV1 rgpW1] = ebcmprod(rgc1.',rgp1.',rgb1.',U1.',V1.',W1.');
    [rgcU2 rgcV2 rgbU2 rgbV2 rgpW2] = ebcmprod(rgc2,rgp2,rgb2,U2,V2,W2);
            
    % Outer products for Q: (l,m) along rows (index 1), (p,q) along columns (index 2)
    M1hat_cross_rgM2 =  n_dot_r*(cV1*rgcU2 - cU1*rgcV2);
    N1hat_cross_rgN2 =  n_dot_r*(bV1*rgbU2 - bU1*rgbV2) ...
                          + n_dot_th*(bU1*rgpW2 - pW1*rgbU2) + n_dot_ph*(pW1*rgbV2 - bV1*rgpW2);
    M1hat_cross_rgN2 =  n_dot_r*(cU1*rgbU2 + cV1*rgbV2) - n_dot_th*(cV1*rgpW2) - n_dot_ph*(cU1*rgpW2);
    N1hat_cross_rgM2 = -n_dot_r*(bU1*rgcU2 + bV1*rgcV2) + n_dot_th*(pW1*rgcV2) + n_dot_ph*(pW1*rgcU2);
    
    % Q 
    Qmm = Qmm + const*M1hat_cross_rgN2 + N1hat_cross_rgM2; 
    Qmn = Qmn + const*M1hat_cross_rgM2 + N1hat_cross_rgN2;
    Qnm = Qnm + const*N1hat_cross_rgN2 + M1hat_cross_rgM2;
    Qnn = Qnn + const*N1hat_cross_rgM2 + M1hat_cross_rgN2;
    
    % Outer products for RgQ: (l,m) along rows (index 1), (p,q) along columns (index 2)
    rgM1hat_cross_rgM2 =  n_dot_r*(rgcV1*rgcU2 - rgcU1*rgcV2);
    rgN1hat_cross_rgN2 =  n_dot_r*(rgbV1*rgbU2 - rgbU1*rgbV2) ...
                            + n_dot_th*(rgbU1*rgpW2 - rgpW1*rgbU2) + n_dot_ph*(rgpW1*rgbV2 - rgbV1*rgpW2);
    rgM1hat_cross_rgN2 =  n_dot_r*(rgcU1*rgbU2 + rgcV1*rgbV2) - n_dot_th*(rgcV1*rgpW2) - n_dot_ph*(rgcU1*rgpW2);
    rgN1hat_cross_rgM2 = -n_dot_r*(rgbU1*rgcU2 + rgbV1*rgcV2) + n_dot_th*(rgpW1*rgcV2) + n_dot_ph*(rgpW1*rgcU2);
    
    % RgQ
    RgQmm = RgQmm + const*rgM1hat_cross_rgN2 + rgN1hat_cross_rgM2;
    RgQmn = RgQmn + const*rgM1hat_cross_rgM2 + rgN1hat_cross_rgN2;
    RgQnm = RgQnm + const*rgN1hat_cross_rgN2 + rgM1hat_cross_rgM2;
    RgQnn = RgQnn + const*rgN1hat_cross_rgM2 + rgM1hat_cross_rgN2;
    
end

% load Q and RgQ
Q(indM,indM) = Qmm;
Q(indM,indN) = Qmn;
Q(indN,indM) = Qnm;
Q(indN,indN) = Qnn;
RgQ(indM,indM) = RgQmm;
RgQ(indM,indN) = RgQmn;
RgQ(indN,indM) = RgQnm; 
RgQ(indN,indN) = RgQnn; 
    
% right matrix divide for matrix inverse to compute T-matrix
T = -RgQ/Q;

% separate block T-matrices
Tmm = T(indM,indM);
Tmn = T(indM,indN);
Tnm = T(indN,indM);
Tnn = T(indN,indN);
end

function [cj pj bj] = ebcmbessel(tot,L,l,l2,k,r,rgstr)
% ebcm helper function to compute Bessel functions
%
% tot:      total number of harmonics
% L:        maximum degree harmonic
% l:        degree l mapped to harmonic linear index (1,1,1,2,2,2,2,2,etc)
% l2:       array of l index 1:L
% k:        wavenumber
% r:        [1x1] radial point
% rgstr:    'rg' for regular form of Bessel functions
%
% cj,pj,bj: [1xtot] radial Bessel functions associated 
%           with C,P,B vector spherical harmonics
%
% Dependencies: sbesselj,sbesselj2,sbesseljp,sbesselh,sbesselhp

kr = k*r;
cj = zeros(1,tot);
pj = zeros(1,tot);
bj = zeros(1,tot);
cjtmp = zeros(1,L);
pjtmp = zeros(1,L);
bjtmp = zeros(1,L);

% evalute bessel functions
if strcmp(rgstr,'rg')
    cjtmp(l2) = sbesselj(l2,kr);        % j_l(kr), handles lim r->0 
    pjtmp(l2) = sbesselj2(l2,kr);       % j_l(kr)/kr, "
    bjtmp(l2) = sbesseljp(l2,kr);       % j'_l(kr), " 
else
    cjtmp(l2) = sbesselh(l2,kr);        % h_l^(1)(kr)
    pjtmp(l2) = sbesselh(l2,kr)./kr;    % h_l^(1)(kr)/kr
    bjtmp(l2) = sbesselhp(l2,kr);       % h'_l^(1)(kr)
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

function [cU cV bU bV pW] = ebcmprod(c,p,b,U,V,W);
% ebcm helper function to compute wave function components
cU = c.*U;
cV = c.*V;
bU = b.*U;
bV = b.*V;
pW = p.*W;
end