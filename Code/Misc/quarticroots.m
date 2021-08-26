function [x1 x2 x3 x4] = quarticroots(a,b,c,d,e)
% Vectorized roots of a quartic polynomial
%
%     ax^4 + bx^3 + cx^2 + dx + e = 0
%
% a,b,c,d,c     Coefficients, any size
%
% x1,x2,x3,x4   Polynomial roots, same size as a,b,c,d,e
%
% Dependencies: cubicroots

% put quartic in standard form
a3 = b./a;
a2 = c./a;
a1 = d./a;
a0 = e./a;

% coefficients of the resolvent cubic
cube_a = ones(size(a3));
cube_b = -a2;
cube_c = a1.*a3 - 4*a0;
cube_d = 4*a2.*a0 - a1.^2 - a3.^2.*a0;

% one real root of the cubic
y1 = cubicroots(cube_a,cube_b,cube_c,cube_d);

% compute the parameters
R = sqrt((1/4)*a3.^2 - a2 + y1);
w = -(1/4)*a3;
T = (3/4)*a3.^2 - 2*a2;
U = (1/4)*(4*a3.*a2-8*a1-a3.^3)./R;
V = 2*sqrt(y1.^2-4*a0);

D = zeros(size(a3));
E = zeros(size(a3));
% condition for R ~= 0
I1 = (R~=0);
D(I1) = sqrt(T(I1)-R(I1).^2+U(I1));
E(I1) = sqrt(T(I1)-R(I1).^2-U(I1));

% condition for R = 0
I2 = ~I1;
D(I2) = sqrt(T(I2)+V(I2));
E(I2) = sqrt(T(I2)-V(I2));

% four roots
x1 = w+R/2+D/2;
x2 = w+R/2-D/2;
x3 = w-R/2+E/2;
x4 = w-R/2-E/2;

