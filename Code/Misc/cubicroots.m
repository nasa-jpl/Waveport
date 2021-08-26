function [x1 x2 x3] = cubicroots(a,b,c,d)
% Vectorized roots of a cubic polynomial
%
%     ax^3 + bx^2 + cx + d = 0
%
% a,b,c,d       Coefficients, any size
%
% x1,[x2,x3]    Polynomial roots, same size as a,b,c,d
%               If nargout==1, function returns one real root, x1

% Q, R, and D parameters
Q = (3*a.*c-b.^2)./(9*a.^2);
R = (9*a.*b.*c-27*a.^2.*d-2*b.^3)./(54*a.^3);
D = Q.^3 + R.^2;
w = -b./(3*a);
S = zeros(size(D));
T = zeros(size(D));

% condition for one real root, and two complex
I1 = (D>=0); 
S(I1) = nthroot(R(I1)+sqrt(D(I1)),3);
T(I1) = nthroot(R(I1)-sqrt(D(I1)),3);

% condition for three real roots
I2 = ~I1; 
rho = nthroot(-Q(I2).^3,2);
theta = acos(R(I2)./rho);
S(I2) = nthroot(rho,3).*exp(1i*theta/3);
T(I2) = nthroot(rho,3).*exp(-1i*theta/3);

% roots 
x1 = real(w + S + T);
x2 = w - (1/2)*(S+T) + (1i*sqrt(3)/2)*(S-T);
x3 = w - (1/2)*(S+T) - (1i*sqrt(3)/2)*(S-T);

