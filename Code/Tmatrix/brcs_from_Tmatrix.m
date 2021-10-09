function [sig_tt sig_tp sig_pt sig_pp] = brcs_from_Tmatrix(Tmm,Tmn,Tnm,Tnn,L,k,theta_s,phi_s,beta_i,beta_s)
% Backscatter radar cross section from T-matrix 
%
% Tmm,Tmn,Tnm,Tnn:      T-matrix block matrices, size [NxN], N = L^2 + 2L
% L:                    Maximum degree harmonic L 
% k:                    Background wavenumber
% theta_s,phi_s:        Back scatttered directions, radians, any size
% beta_i,beta_s:        [optional] Polarization angles, radians, [1 x 1]
%
% sig_tt,sig_tp         Backscatter radar cross section
% sig_pt,sig_pp:        - Use nargout == 1 for an arbitrary polarization and then 
%                           output is in sig_tt and inputs beta_i and beta_s are needed
%                       - Use nargout == 4 for all four orthogonal polarization combiations and then 
%                           inputs beta_i and beta_s are not needed
%                     
% Dependencies:         BC, lmtable

% size of scattered direction array
sz_s = size(theta_s);

% total number of scattered directions
Ns = numel(theta_s);

% reshape the inputs
theta_s = theta_s(:);
phi_s = phi_s(:);

% apply coefficient matricies as sparse diagonal matricies to the T-matrix
N = L^2 + 2*L;
tab = lmtable(L);
l = tab(:,1);
l1 = (-1i).^(l+1);
l2 = (-1i).^(l);
diagind = 1:N;
L1 = sparse(diagind,diagind,l1,N,N);
L2 = sparse(diagind,diagind,l2,N,N);
Tmm = L1*Tmm*L2;
Tmn = L1*Tmn*(-L1);
Tnm = L2*Tnm*L2;
Tnn = L2*Tnn*(-L1);

% scattered direction vector spherical harmonics
[Bth_s, Bphi_s, Cth_s, Cphi_s] = BC(L,theta_s,phi_s,'norm');

% build the matrices and multiply (essentially the S-matrix)
S = zeros(Ns,2,2);
for n=1:Ns,
    S(n,:,:) = ([Cth_s(n,:) Bth_s(n,:); Cphi_s(n,:) Bphi_s(n,:)]* ...
        ([Tmm Tmn; Tnm Tnn]*[Cth_s(n,:)' Cphi_s(n,:)'; Bth_s(n,:)' Bphi_s(n,:)']));
end

% Extract the four components and reshape
Stt = reshape(S(:,1,1),sz_s);
Stp = reshape(S(:,1,2),sz_s);
Spt = reshape(S(:,2,1),sz_s);
Spp = reshape(S(:,2,2),sz_s);

const = (4*pi)^3/abs(k)^2;
if nargout == 1 % arbitrary polarization
   sig_tt = const*abs(cos(beta_s)*Stt*cos(beta_i) + cos(beta_s)*Stp*sin(beta_i) ...
       + sin(beta_s)*Spt*cos(beta_i) + sin(beta_s)*Spp*sin(beta_i)).^2;   
elseif nargout == 4 % four component polarizations
   sig_tt = const*abs(Stt).^2;
   sig_tp = const*abs(Stp).^2;
   sig_pt = const*abs(Spt).^2;
   sig_pp = const*abs(Spp).^2;    
else
    error('bad output')
end