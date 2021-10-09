function [Stt, Stp, Spt, Spp] = compute_S_from_T(Tmm,Tmn,Tnm,Tnn,L,k,theta_s,phi_s,theta_i,phi_i)
% Compute S-matrix points given full T-matrix 
%
% Tmm,Tmn,Tnm,Tnn:    T-matrix block matrices
%                     Size [NxN], N = L^2 + 2L
% L:                  Maximum degree harmonic L 
% k:                  Background wavenumber
% theta_s,phi_s:      Scatttered directions, radians, any size, sz_s = size(theta_s)
% theta_i,phi_i:      Incident directions, radians, any size, sz_i = size(theta_i)
%
% Stt,Stp,Spt,Spp:    S-matrix block matrices Theta/Phi components, 
%                     Size [sz_s sz_i], concatenated sizes of the input
%                     incident and scattered directions
%                     
% Dependencies:       BC, lmtable

% size of incident direction arrays
sz_i = size(theta_i);
% size of scattered direction arrays
sz_s = size(theta_s);

% total number of incident/scattered directions
Ni = numel(theta_i);
Ns = numel(theta_s);

% reshape the inputs
theta_i = theta_i(:);
phi_i = phi_i(:);
theta_s = theta_s(:);
phi_s = phi_s(:);

% apply coefficient matricies as sparse diagonal matricies to the T-matrix
N = L^2 + 2*L;
tab = lmtable(L);
l = tab(:,1);
l1 = (1i).^(-l-1);
l2 = (1i).^(-l);
diagind = 1:N;
L1 = sparse(diagind,diagind,l1,N,N);
L2 = sparse(diagind,diagind,l2,N,N);
L1inv = -sparse(diagind,diagind,1./l1,N,N);
L2inv = sparse(diagind,diagind,1./l2,N,N);
Tmm = L1*Tmm*L2inv;
Tmn = L1*Tmn*L1inv;
Tnm = L2*Tnm*L2inv;
Tnn = L2*Tnn*L1inv;

% incident direction vector spherical harmonics
[Bth_i, Bphi_i, Cth_i, Cphi_i] = BC(L,theta_i,phi_i,'norm');

% scattered direction vector spherical harmonics
[Bth_s, Bphi_s, Cth_s, Cphi_s] = BC(L,theta_s,phi_s,'norm');

% build the matrices and multiply
Smatrix = (4*pi/k)*([Cth_s Bth_s; Cphi_s Bphi_s]*...
    ([Tmm Tmn; Tnm Tnn]*[Cth_i' Cphi_i'; Bth_i' Bphi_i']));

% indecies for different components in the full matrix
grab_i_theta = 1:Ni;
grab_s_theta = 1:Ns;
grab_i_phi = (Ni+1):(2*Ni);
grab_s_phi = (Ns+1):(2*Ns);

% Extract the four components and reshape
Stt = reshape(Smatrix(grab_s_theta,grab_i_theta),[sz_s sz_i]);
Stp = reshape(Smatrix(grab_s_theta,grab_i_phi),[sz_s sz_i]);
Spt = reshape(Smatrix(grab_s_phi,grab_i_theta),[sz_s sz_i]);
Spp = reshape(Smatrix(grab_s_phi,grab_i_phi),[sz_s sz_i]);

