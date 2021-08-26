function [sn] = SN(N,theta,thetao)
% Windowing function for FMM translation operator interpolation
%
% S_N(theta,theta_o) = R_N(theta,theta_o)/R_N(0,theta_o) 
%
% R_N(theta,theta_o) = sinh((2N+1)sinh^-1 sqrt(sin^2(theta_o/2)-sin^2(theta/2)))/...
%                           sqrt(sin^2(theta_o/2)-sin^2(theta/2))
%
% N:        window width parameter
% theta:    interpolation point, angle in radians
% thetao: 	stepping parameter 
%
% sn:       window function

sqtheta = sqrt(sin(thetao/2).^2 - sin(theta/2).^2);
RNtheta = sinh((2*N+1)*asinh(sqtheta))./sqtheta;
sqo = sin(thetao/2);
RNo = sinh((2*N+1).*asinh(sqo))./sqo;
sn = RNtheta./RNo;
ind = find(theta == -thetao);
if ~isempty(ind)
    sn(ind) = 1;
end