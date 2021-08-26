function [RX RY RZ theta1 theta2 l1 l2] = refractionPlane(x1,y1,z1,X2,Y2,Z2,zo,n1,n2,str);
% Refraction through a dieletric interface
% Interface is parallel to XY plane. Assumes n1, n2 are real.
% Default solution is via quartic, with optional iterative solution
%
% x1,y1,z1  1x1 Cartesian coordinates of point in medium 1
% X2,Y2,Z2  [any size] Cartesian coordinates of points in medium 2
% zo        z coordinate of dielectric interface
% n1,n2     Index of refraction in medium 1 and 2
% str       [optional] 'it' for iterative solution
%
% RX,RY,RZ  [size(X2)] Cartesian coordinates of refraction point
% theta1    [size(X2)] Refraction angle in medium 1 (radians)
% theta2    [size(X2)] Refraction angle in medium 2 (radians)
% l1        [size(X2)] Distance from r1 to r
% l2        [size(X2)] Distance from r2 to r
%
% Dependencies: iterateRefractionTest, iterateRefraction

sz = size(X2);
N = numel(X2);
n1 = real(n1);
n2 = real(n2);

if n1 < 1 || n2 < 1
    disp('Error: bad index of refraction')
    return
end
if z1 <= zo  
    disp('Error: bad z coordinate')
    return
end

% projection onto vertical plane
x21 = X2 - x1;
y21 = Y2 - y1;
L = sqrt(x21.^2 + y21.^2);
uhat_x = x21./L;
uhat_y = y21./L;
H = abs(z1 - zo);
D = abs(Z2 - zo);

if nargin == 9              % quartic solution
    % quartic coefficients 
    a = (n1.^2 - n2.^2)*ones(sz);
    b = -2*L.*a;
    c = n1.^2.*(L.^2 + D.^2) - n2.^2.*(L.^2 + H.^2);
    d = 2*n2.^2.*L.*H.^2;
    e = -n2.^2.*L.^2.*H.^2;
    % quartic roots vectorized
    [w1 w2 w3 w4] = quarticroots(a,b,c,d,e);
    % pick out the one real, positive, valid root
    r = [w1(:) w2(:) w3(:) w4(:)];
    u = zeros(N,1);
    for n=1:4,
        ind = and(and(imag(r(:,n))==0,r(:,n)>0),r(:,n)<L(:));
        u(ind) = r(ind,n);
    end
    u = reshape(u,sz);
elseif nargin == 10 && strcmp(str,'it')   % iterative solution 
    % maximum straight ray theta1
    [M I] = max(atan2(L,abs(Z2-z1)));
    % test for number of iterations
    nit = iterateRefractionTest(L(I),H,D(I),n1,n2);
    % use nit iterations for all points
    u = iterateRefraction(L,H,D,n1,n2,nit);
end

% reproject the refraction point in original frame
RX = x1 + uhat_x.*u;
RY = y1 + uhat_y.*u;
RZ = zo*ones(sz);

% compute angles and lengths
theta1 = real(asin(u./sqrt(H.^2 + u.^2)));
theta2 = real(asin((L - u)./(D.^2 + (L - u).^2)));
l1 = sqrt(u.^2 + H.^2);
l2 = sqrt(D.^2 + (L - u).^2);

% if subsurface points are vertically algined with upper point
% expect this number to be small so use find
% other parameters are taken care of by L
ind = find(L==0);
RX(ind) = x1;
RY(ind) = y1;
theta1(ind) = 0;
theta2(ind) = 0;
l1(ind) = H;
l2(ind) = D(ind);

% take care of points above or equal to boundary
ind = find(Z2 >= zo);
RX(ind) = nan;
RY(ind) = nan;
RZ(ind) = nan;
theta1(ind) = real(asin(L(ind)./sqrt((H-Z2(ind)).^2 + L(ind).^2)));
theta2(ind) = nan;
l1(ind) = sqrt(L(ind).^2 + (H-Z2(ind)).^2);
l2(ind) = nan;

% take care of points on boundary
ind = find(Z2 == zo);
RX(ind) = X2(ind);
RY(ind) = Y2(ind);
RZ(ind) = Z2(ind);
end


function [it] = iterateRefractionTest(L,h,d,n1,n2);
% Testing function that returns the number of iterations for one point
    th1 = atan2(L,h+d);
    it = 0;
    ratio = 1;
    while ratio > 1e-8
        th2 = asin((n1/n2)*sin(th1));
        u = L-d.*tan(th2);
        th1new = atan2(u,h);        
        ratio = abs(th1new-th1);
        th1 = th1new;
        it = it+1;
    end
end

function [u] = iterateRefraction(L,h,d,n1,n2,nit);
% Vectorized iteration solution for given number of iterations, nit
    th1 = atan2(L,h+d);
    u = d.*tan(th1);
    for n=1:nit
        th2 = asin((n1/n2)*sin(th1));
        u = L-d.*tan(th2);
        th1 = atan2(u,h);
    end
end

