function [Tmm, Tmn, Tnm, Tnn] = convert_S_to_T(Stt,Stp,Spt,Spp,k,I,J,L)
% S-matrix to T-matrix transformation using quadrature
%
% Stt,Stp,Spt,Spp:    S-matrix block matrices Theta/Phi components, 
%                     Size [IxJxIxJ], sampled on the points of quadrature
%                     Scattered directions dimensions 1 and 2
%                     Incident directions dimensions 3 and 4
% k:                  Background wavenumber
% I,J:                I = 2L+1, J = L+1,
% L:                  Maximum degree harmonic L
%     
% Tmm,Tmn,Tnm,Tnn:    T-matrix block matrices
%                     Size [NxN], N = L^2 + 2L
%                     
% Dependencies:       vst, lmtable, legpts, Plm, Plmp

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
S2 = zeros(2*Nk,2*N);
Tmatrix = zeros(2*N);

%%% compute the right-hand multiplication over incident field directions
% (rows of the S-matrix) using the conjugate transpose for forward
% spherical transform

% for each scattered direction, extract the incident directions as a field
for n=1:2*Nk,
    % top half of full S-matrix
    if n <= Nk
        [indi indj] = ind2sub([I J],n);
        Fth = squeeze(Stt(indi,indj,:,:));
        Fphi = squeeze(Stp(indi,indj,:,:));
    else % bottom half of full S-matrix
        [indi indj] = ind2sub([I J],n-Nk);
        Fth = squeeze(Spt(indi,indj,:,:));
        Fphi = squeeze(Spp(indi,indj,:,:));    
    end
    % compute forward transform on conjugate field
    [blm clm] = vst(conj(Fth),conj(Fphi),L,muj,wj,'norm',plm,dplm);
     
    % conjgate transpose to bring it back and store in temporary array
    S2(n,grab1_lm) = conj(clm.');
    S2(n,grab2_lm) = conj(blm.');   
end

%%% compute the left-hand multiplication over the scattered directions 
% this indexes over the harmonics of the first transform (rows of the
% temporary array)
for n=1:2*N,
    Fth = reshape(S2(grab1_k,n),I,J);
    Fphi = reshape(S2(grab2_k,n),I,J);
    
    % apply forward transform
	[blm clm] = vst(Fth,Fphi,L,muj,wj,'norm',plm,dplm);

    % put the transform down the columns of the Tmatrix
    Tmatrix(grab1_lm,n) = clm;
    Tmatrix(grab2_lm,n) = blm;   
end

% build coefficient matricies as sparse diagonal matricies
tab = lmtable(L);
l = tab(:,1);
l1 = (1i).^(-l-1);
l2 = (1i).^(-l);
diagind = 1:(2*N);
L1invL2inv = sparse(diagind,diagind,1./[l1; l2],2*N,2*N);
L2L1 = sparse(diagind,diagind,[l2; -l1],2*N,2*N);

% apply coefficient matricies
Tmatrix = (k/(4*pi))*(L1invL2inv*(Tmatrix*L2L1));

% Extract the four components
Tmm = Tmatrix(grab1_lm,grab1_lm);
Tmn = Tmatrix(grab1_lm,grab2_lm);
Tnm = Tmatrix(grab2_lm,grab1_lm);
Tnn = Tmatrix(grab2_lm,grab2_lm);

