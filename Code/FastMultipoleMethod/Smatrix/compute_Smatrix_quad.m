function [Fth, Fphi] = compute_Smatrix_quad(Gth,Gphi,Stt,Stp,Spt,Spp,I,J)
% Application of S-matrix to field pattern using quadrature
%
% Gth,Gphi:           [IxJ] incident field pattern sampled at the nodes of quadrature
% Stt,Stp,Spt,Spp:    S-matrix block matrices Theta/Phi components, 
%                     Size [IxJxIxJ], sampled on the points of quadrature
%                     Scattered directions dimensions 1 and 2
%                     Incident directions dimensions 3 and 4
% I,J:                I = 2L+1, J = L+1,
%     
% Fth,Fphi:           [IxJ] scattered field pattern sampled at nodes of quadrature
%                     
% Dependencies:       legpts

% compute quadrature weights
[muj, wj] = legpts(J);

% apply quadrature weights to incident field pattern
const = 2*pi/I;
W = repmat(wj(:)',I,1);
Gth = const*W.*Gth;
Gphi = const*W.*Gphi;

% matrix multiplication
N = I*J;
Fth = reshape(Stt,N,N)*Gth(:) + reshape(Stp,N,N)*Gphi(:);
Fphi = reshape(Spt,N,N)*Gth(:) + reshape(Spp,N,N)*Gphi(:);
Fth = reshape(Fth,I,J);
Fphi = reshape(Fphi,I,J);

