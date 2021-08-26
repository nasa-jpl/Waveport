function [alpha beta gamma] = rot2euler(R)
% Convert a 3x3 rotation matrix to the Euler angles alpha, beta, gamma
%
% R:                3x3 rotation matrix
%
% alpha,beta,gamma: ZXZ Euler angles (randians)

    alpha = atan2(R(1,3),-R(2,3));
    beta = atan2(sqrt(R(1,3)^2+R(2,3)^2),R(3,3));
    gamma = atan2(R(3,1),R(3,2));
end