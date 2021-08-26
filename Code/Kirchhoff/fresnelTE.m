function [RTE] = fresnelTE(er1,er2,theta_inc_deg)
% Fresnel reflection coefficient - TE
%
% er1, er2:         Relative permittivity (complex) in regions 1 and 2
% theta_inc_deg:    Incident angle (degrees)
%
% RTE:              TE reflection coefficient

ratio = er2./er1;
sq = sqrt(ratio-sind(theta_inc_deg).^2);
RTE = (cosd(theta_inc_deg) - sq)./(cosd(theta_inc_deg) + sq);

