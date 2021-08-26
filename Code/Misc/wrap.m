function [Phi] = wrap(phi)
% Phase wrapping with analytical functions
%
% phi:  phase (rad)
%
% Phi:  wrapped phase (rad) on [-pi pi]

    Phi = -2*atan(cot(phi/2 + pi/2));
end