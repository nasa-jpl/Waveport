function [Tbb, Tbc, Tcb, Tcc] = convert_T_to_Tfar(Tmm,Tmn,Tnm,Tnn,L,k)
% Convert T-matrix to far-field T-matrix
%
% Tmm,Tmn,Tnm,Tnn:    T-matrix block matrices
%                     Size [NxN], N = L^2 + 2L
% L:                  Maximum degree harmonic L 
% k:                  Background wavenumber
%
% Tbb,Tbc,Tcb,Tcc:    T-matrix block matrices
%                     Size [NxN], N = L^2 + 2L
%                     
% Dependencies:       lmtable

% apply coefficient matrices as sparse diagonal matrices to the T-matrix
N = L^2 + 2*L;
tab = lmtable(L);
l = tab(:,1);
l1 = (1i).^(-l-1);
l2 = (1i).^(-l);
diagind = 1:N;
L1 = sparse(diagind,diagind,l1,N,N);
L2 = sparse(diagind,diagind,l2,N,N);
L1inv = sparse(diagind,diagind,1./l1,N,N);
L2inv = sparse(diagind,diagind,1./l2,N,N);
const = 4*pi/k;
Tbb = (const*L2)*Tnn*L1inv;
Tbc = (const*L2)*Tnm*L2inv;
Tcb = (const*L1)*Tmn*L1inv;
Tcc = (const*L1)*Tmm*L2inv;

