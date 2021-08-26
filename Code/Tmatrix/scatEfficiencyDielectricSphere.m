function [Qsca] = scatEfficiencyDieletricSphere(k1,k2,a)
% Scattering Efficiency of a Dieleitric Sphere
%
% k1:   background wave number
% k2:   sphere wave number   
% a:    sphere radius
%
% Dependencies: scatCrossSectionDieletricSphere 

[sig_sca] = scatCrossSectionDieletricSphere(k1,k2,a);
sig_geo = pi*a.^2;
Qsca = sig_sca./sig_geo;







