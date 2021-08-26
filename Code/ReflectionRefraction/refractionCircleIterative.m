function [rx ry theta1 theta2 thetal v1mag v2mag] = refractionCircleIterative(x1,y1,x2,y2,n1,n2,r,nit)
% Refraction through a circular interface, iterative solution
% assumes n1, n2 are real
%
% x1,y1     Coordinates of exterior point, r1, centered on circle
% x2,y2     Coordinates of interior point, r2, centered on circle
% n1        Index of refraction for exterior
% n2        Index of refraction for interior
% r         Radius circle
% nit       Number of iterations
%
% rx,ry         Coordinates of refraction point
% theta1        Incidence angle
% theta2        Transmission angle
% thetal        Look angle measured from perpendicular
% v1mag, v2mag  Lengths of the vectors from each point to the refraction point

sz =size(x1);
n1 = real(n1);
n2 = real(n2);

% rotate points to align r_1 with y-axis
x1 = x1(:);
y1 = y1(:);
N = length(x1);
for n=1:N,
    theta_rot = atan2(y1(n),x1(n));
    theta_rot = -(theta_rot - pi/2);
    R = [cos(theta_rot) -sin(theta_rot); 
       sin(theta_rot) cos(theta_rot)];
    r1rot = R*[x1(n)'; y1(n)'];
    x1(n) = r1rot(1);
    y1(n) = r1rot(2);
    r2rot = R*[x2(n); y2(n)];
    x2(n) = r2rot(1);
    y2(n) = r2rot(2);
end

% aux parameters
h = y1 - r;
d = r - y2;

% default number o iterations
if nargin == 7
    nit = 40;
end

% initialized iterations
thetal = atan2(x2,h + d);

% find points along radial line
ind = find(thetal==0);

for n=1:nit
    % compute incident, transmission and point angles
    sin_theta1 = asin(((r + h)./r).*sin(thetal));
    theta1 = asin(sin_theta1);
    theta2 = asin(n1./n2.*sin_theta1);
    theta_p = theta2 - (theta1 - thetal);

    % solve for intersection point
    m = -cot(theta_p);
    b = x2.*cot(theta_p) + y2;
    A = m.^2 + 1;
    B = 2*m.*b;
    C = -r.^2 + b.^2;
    rx = (-B-sign(x2).*sqrt(B.^2 - 4*A.*C))./(2*A);
    ry = (m.*rx + b);

    % if refraction point is on wrong side of sphere, use other solution
    ind2 = (ry<y2);
    rx(ind2) = (-B(ind2) + sign(x2(ind2)).*sqrt(B(ind2).^2 - 4*A(ind2).*C(ind2)))./(2*A(ind2));
    ry(ind2) = (m(ind2).*rx(ind2) + b(ind2));    
    
    % update look angle
    thetal = atan2(rx,h + r - ry);
end

% take care of radial points
rx(ind) = 0;
ry(ind) = r;
thetal(ind) = 0;
theta1(ind) = 0;
theta2(ind) = 0;



% compute vectors for remaining solutions
v1 = [(x1-rx), (y1-ry)];
v2 = [(x2-rx), (y2-ry)];
v1mag = sqrt(sum(v1.^2,2));
v2mag = sqrt(sum(v2.^2,2));

% rotate solution back
for n=1:N,
    theta_rot = atan2(y1(n),x1(n));
    theta_rot = -(theta_rot - pi/2);
    R = [cos(theta_rot) -sin(theta_rot); 
       sin(theta_rot) cos(theta_rot)];
    rrrot = R'*[rx(n); ry(n)];
    rx(n) = rrrot(1);
    ry(n) = rrrot(2);
end

rx=reshape(rx,sz);
ry=reshape(ry,sz);
theta1=reshape(theta1,sz);
theta2=reshape(theta2,sz);
thetal=reshape(thetal,sz);
v1mag=reshape(v1mag,sz);
v2mag=reshape(v2mag,sz);



