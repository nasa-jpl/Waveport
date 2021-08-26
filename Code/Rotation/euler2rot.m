function [R] = euler2rot(alpha,beta,gamma)
% 3D rotation matrix for ZXZ Euler angles
%
% alpha:    right-handed rotation angle about Z (radians)
% beta:     right-handed rotation angle about X (radians)
% gamma:    right-handed rotation angle about Z (radians)
%
% R:        [3x3] rotation matrix

    Rza = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
    Rxb = [1 0 0; 0 cos(beta) -sin(beta); 0 sin(beta) cos(beta)];
    Rzg = [cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];	
    R = Rza*Rxb*Rzg;
end
