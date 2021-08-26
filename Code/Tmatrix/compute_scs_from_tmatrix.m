function [scs] = compute_scs_from_tmatrix(Tmm,Tmn,Tnm,Tnn,L,k,theta_i,phi_i,beta)
% Compute scattering cross section from T-matrix
%
% Tmm,Tmn,Tnm,Tnn:    T-matrix block matrices
%                     Size [NxN], N = L^2 + 2L
% L:                  Maximum degree harmonic L 
% k:                  Background wavenumber
% theta_i,phi_i:      Incident directions, radians, any size, sz_i = size(theta_i)
% beta:               Right-hand polarization angle (radians) relative to
%                     \hat{\theta}, any size, sz_b = size(beta)
%
% scs:                Scattering cross section 
%                     Size [sz_i sz_b], concatenated sizes of the incident
%                     directions and beta (singular dimensions squeezed)
%                     
% Dependencies:       BC, lmtable

% size of incident direction arrays and incident directions
sz_i = size(theta_i);
sz_b = size(beta);

% total number of incident directions
Ni = numel(theta_i);
Nb = numel(beta);
scs = zeros(Ni,Nb);

% reshape the inputs
theta_i = theta_i(:);
phi_i = phi_i(:);
beta = beta(:);

% incident direction vector spherical harmonics (sized [Ni x N])
[Bth_i, Bphi_i, Cth_i, Cphi_i] = BC(L,theta_i,phi_i,'norm');

% create coefficient matricies as sparse diagonal matricies 
N = L^2 + 2*L;
tab = lmtable(L);
l = tab(:,1);
l1 = (1i).^(l);
l2 = (1i).^(l+1);
diagind = 1:N;
L1 = sparse(diagind,diagind,l1,N,N);
L2 = sparse(diagind,diagind,l2,N,N);

% separate the polarizations, apply conjugate, perform sum of l'm' with
% matrix multiplication. the results are sized [N x Ni].
Th1 = (Tmm*L1)*Cth_i' + (Tmn*L2)*Bth_i';
Ph1 = (Tmm*L1)*Cphi_i' + (Tmn*L2)*Bphi_i';
Th2 = (Tnm*L1)*Cth_i' + (Tnn*L2)*Bth_i';
Ph2 = (Tnm*L1)*Cphi_i' + (Tnn*L2)*Bphi_i';

% for each \beta angle, compute abs^2 and sum over lm
for b=1:Nb,
    tmp1 = sum(abs(Th1*cos(beta(b)) + Ph1*sin(beta(b))).^2,1);
    tmp2 = sum(abs(Th2*cos(beta(b)) + Ph2*sin(beta(b))).^2,1);
    scs(:,b) = tmp1.' + tmp2.';
end

% apply constant 
scs = (16*pi^2/(abs(k)^2))*scs;

% reshape and squeeze singular dimensions
scs = squeeze(reshape(scs,[sz_i sz_b]));



