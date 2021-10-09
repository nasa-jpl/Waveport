function [Stt, Stp, Spt, Spp] = convert_T_to_S(Tmm,Tmn,Tnm,Tnn,k,I,J,L)
% T-matrix to S-matrix transformation using quadrature
%
% Tmm,Tmn,Tnm,Tnn:    T-matrix block matrices
%                     Size [NxN], N = L^2 + 2L
% k:                  Background wavenumber
% I,J:                I = 2L+1, J = L+1, S-matrix quadratute sampling
% L:                  Maximum degree harmonic L 
%
% Stt,Stp,Spt,Spp:    S-matrix block matrices Theta/Phi components, 
%                     Size [IxJxIxJ], sampled on the points of quadrature
%                     Scattered directions dimensions 1 and 2
%                     Incident directions dimensions 3 and 4
%                     
% Dependencies:       ivst, lmtable, legpts, Plm, Plmp

% precompute the Legendre polynomials for the transforms
[muj, wj] = legpts(J);
plm = Plm(L,muj); % includes monopole
dplm = Plmp(L,muj); % includes monopole

% indexing
Nk = I*J;
grab1_k = 1:Nk;
grab2_k = (Nk+1):(2*Nk);
N = L^2 + 2*L;
grab1_lm = 1:N;
grab2_lm = (N+1):(2*N);

% temporary arrays for storing intermediate transform and full T-matrix
S2 = zeros(2*N,2*Nk);
Smatrix = zeros(2*Nk,2*Nk);

%%% compute the right-hand multiplication over incident harmonics
% (rows of the T-matrix) using the conjugate transpose for the inverse
% spherical transform for each harmonic, extract the scattered harmonics
% apply normalization if 

% apply coefficients as sparse diagonal matricies
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

for n=1:2*N,
    % top half of full T-matrix
    if n <= N
        clm = squeeze(Tmm(n,:));
        blm = squeeze(Tmn(n,:));
    else % bottom half of full T-matrix
        clm = squeeze(Tnm(n-N,:));
        blm = squeeze(Tnn(n-N,:));
    end
    % compute inverse transform on conjugate coefficients
    [Fth Fphi] = ivst(conj(blm(:)),conj(clm(:)),L,muj,'norm',plm,dplm);
    % conjugate transpose to bring it back and store in temporary array
    S2(n,grab1_k) = conj(Fth(:).');
    S2(n,grab2_k) = conj(Fphi(:).');   
end

%%% compute the left-hand multiplication over the scattered harmonics
% this indexes over the incident field points of the first transform (rows of the
% temporary array)
for n=1:2*Nk,
    clm = S2(grab1_lm,n);
    blm = S2(grab2_lm,n);
    % compute inverse transform
	[Fth Fphi] = ivst(blm,clm,L,muj,'norm',plm,dplm);
    % put the scattered field values down the columns of Smatrix
    Smatrix(grab1_k,n) = Fth(:);
    Smatrix(grab2_k,n) = Fphi(:);   
end

sz = [I J I J];
% Extract the four components and reshape
Stt = (4*pi/k)*reshape(Smatrix(grab1_k,grab1_k),sz);
Stp = (4*pi/k)*reshape(Smatrix(grab1_k,grab2_k),sz);
Spt = (4*pi/k)*reshape(Smatrix(grab2_k,grab1_k),sz);
Spp = (4*pi/k)*reshape(Smatrix(grab2_k,grab2_k),sz);