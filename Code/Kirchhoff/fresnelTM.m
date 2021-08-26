function [RTM] = fresnelTM(er1,er2,theta_inc_deg)
% Fresnel reflection coefficient - TM
%
% er1, er2:         Relative permittivity (complex) in regions 1 and 2
% theta_inc_deg:    Incident angle (degrees)
%
% RTM:              TM reflection coefficient

ratio = er2./er1;
sq = sqrt(ratio-sind(theta_inc_deg).^2);
RTM = (ratio.*cosd(theta_inc_deg) - sq)./(ratio.*cosd(theta_inc_deg) + sq);
