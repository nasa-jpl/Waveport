function [Svv Svh Shv Shh] = smatrix_thin_circular_disk(k,a,L,er,theta_s,phi_s,theta_i,phi_i)
% S-matrix of a thin circular disk, oriented with normal along the z-axis
%
% k:                  Background wavenumber (1/m)
% a:                  disk radius (m)
% L:                  disk thickness (m)
% er:                 disk complex relative permittivity 
% theta_s,phi_s:      Scatttered directions, radians, any size, sz_s = size(theta_s)
% theta_i,phi_i:      Incident directions, radians, any size, sz_i = size(theta_i)
% 
% Svv,Svh,Shv,Shh:    S-matrix block matrices h/v components, 
%                     Size [sz_s sz_i], concatenated sizes of the input
%                     incident and scattered directions
%
% Dependencies: jinc

% size of incident/scattered direction arrays
sz_i = size(theta_i);
sz_s = size(theta_s);

% total number of incident/scattered directions
Ni = numel(theta_i);
Ns = numel(theta_s);

% reshape the inputs
theta_i = theta_i(:);
phi_i = phi_i(:);
theta_s = theta_s(:);
phi_s = phi_s(:);

% ndgrid reshaped inputs
[Th_s Th_i] = ndgrid(theta_s,theta_i);
[Phi_s Phi_i] = ndgrid(phi_s,phi_i);

% polarization vector
Pxx = er - 1;
Pyy = Pxx;
Pzz = (er - 1)./er;

% C matrix
sth_s = sin(Th_s);
cth_s = cos(Th_s);
sph_s = sin(Phi_s);
cph_s = cos(Phi_s);
sth_i = sin(Th_i);
cth_i = cos(Th_i);
sph_i = sin(Phi_i);
cph_i = cos(Phi_i);
Cvv = Pzz.*sth_s.*sth_i + Pxx.*cth_s.*cph_s.*cth_i.*cph_i + Pyy.*cth_s.*cth_i.*sph_s.*sph_i;
Cvh = cth_s.*(Pyy.*sph_s.*cph_i - Pxx.*cph_s.*sph_i);
Chv = cth_i.*(Pyy.*cph_s.*sph_i - Pxx.*sph_s.*cph_i);
Chh = Pyy.*cph_s.*cph_i + Pxx.*sph_s.*sph_i;

% constant
const = k^2*L*a^2/2;

% sinc function and constant
Omega2 = (sth_i.*cph_i - sth_s.*cph_s).^2 + (sth_i.*sph_i - sth_s.*sph_s).^2;
U = k*L*sqrt(Omega);
jincU = const*jinc(U);  % matlab's sinc is sin(pi*x)/(pi*x);

% Smatrix elements and reshape
Svv = squeeze(reshape(jincU.*Cvv,[sz_s sz_i]));
Svh = squeeze(reshape(jincU.*Cvh,[sz_s sz_i]));
Shv = squeeze(reshape(jincU.*Chv,[sz_s sz_i]));
Shh = squeeze(reshape(jincU.*Chh,[sz_s sz_i]));




