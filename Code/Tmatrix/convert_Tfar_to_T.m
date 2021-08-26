function [Tmm, Tmn, Tnm, Tnn] = convert_Tfar_to_T(Tbb,Tbc,Tcb,Tcc,L,k)
% Convert far-field T-matrix to T-matrix
%
% Tbb,Tbc,Tcb,Tcc:    T-matrix block matrices
%                     Size [NxN], N = L^2 + 2L
% L:                  Maximum degree harmonic L 
% k:                  Background wavenumber
%
% Tmm,Tmn,Tnm,Tnn:    T-matrix block matrices
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
const = k/(4*pi);
Tnn = (const*L2inv)*Tbb*L1;
Tnm = (const*L2inv)*Tbc*L2;
Tmn = (const*L1inv)*Tcb*L1;
Tmm = (const*L1inv)*Tcc*L2;

