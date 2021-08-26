function [rx ry theta1 theta2 thetal v1mag v2mag] = refractionCircle(x1,y1,x2,y2,n1,n2,r)
% Refraction through a circular interface
% assumes n1, n2 are real
%
% x1,y1     Coordinates of exterior point, r1, centered on circle
% x2,y2     Coordinates of interior point, r2, centered on circle
% n1        Index of refraction for exterior
% n2        Index of refraction for interior
% r         Radius circle
%
% rx,ry         Coordinates of refraction point
% theta1        Incidence angle
% theta2        Transmission angle
% thetal        Look angle measured from perpendicular
% v1mag, v2mag  Lengths of the vectors from each point to the refraction point

if sqrt(x1^2 + y1^2) < r
    error('exterior point is interior')
end
if sqrt(x2^2 + y2^2) > r
    error('interior point is exterior')
end
n1 = real(n1);
n2 = real(n2);
% rotate points to align r_1 with y-axis
theta_rot = atan2(y1,x1);
theta_rot = -(theta_rot - pi/2);
R = [cos(theta_rot) -sin(theta_rot); 
       sin(theta_rot) cos(theta_rot)];
r1rot = R*[x1; y1];
x1 = r1rot(1);
y1 = r1rot(2);
r2rot = R*[x2; y2];
x2 = r2rot(1);
y2 = r2rot(2);

% compute c coefficients
c_1 = (n1^2/n2^2)*y1^2;
c_2 = r^2 + y1^2;
c_3 = -2*y1;
c_4 = y2^2 - x2^2;
c_5 = -2*x2*y2;
c_6 = r^2*x2^2;
c_7 = x2^2 + y2^2 + r^2;
c_8 = -2*x2;
c_9 = -2*y2;

% compute b coefficients
b_1 = c_3*c_5 + c_1*c_8;
b_2 = c_1*c_7 - c_2*c_4;
b_3 = -c_3*c_5*r^2;
b_4 = -c_2*c_6;
b_5 = c_3*c_4 - c_1*c_9;
b_6 = c_2*c_5;
b_7 = c_3*c_6;

% compute a coefficients
a_6 = b_1^2 + b_5^2;
a_5 = 2*b_1*b_2 + 2*b_5*b_6;
a_4 = b_2^2 - b_5^2*r^2 + 2*b_7*b_5 + b_6^2 + 2*b_1*b_3;
a_3 = -2*b_5*b_6*r^2 + 2*b_1*b_4 + 2*b_2*b_3 + 2*b_6*b_7;
a_2 = b_3^2 - b_6^2*r^2 + b_7^2 - 2*b_5*b_7*r^2 + 2*b_2*b_4;
a_1 = -2*b_6*b_7*r^2 + 2*b_3*b_4;
a_0 = b_4^2 - b_7^2*r^2;

% solve for the roots
coef = [a_6 a_5 a_4 a_3 a_2 a_1 a_0];
rt = roots(coef);

% reject roots with imaginary component
ind = (imag(rt)==0);
rx = rt(ind);

% compute ry
ry = sqrt(r^2 - rx.^2);

% compute vectors for remaining solutions
Nr = length(rx);
rr = [rx ry zeros(Nr,1)];
r1 = repmat([0 y1 0],Nr,1);
r2 = repmat([x2 y2 0],Nr,1);
v1 = r1 - rr;
v2 = r2 - rr;
rmag = sqrt(sum(rr.^2,2));
v1mag = sqrt(sum(v1.^2,2));
v2mag = sqrt(sum(v2.^2,2));
sint1 = cross(rr,v1,2)./v1mag./rmag;
sint2 = cross(-rr,v2,2)./v2mag./rmag;
sint1 = sint1(:,3);
sint2 = sint2(:,3);

% compute theta1 and theta2 via dot product
cost1 = dot(rr,v1,2)./v1mag./rmag;
theta1 = acos(cost1);
cost2 = dot(-rr,v2,2)./v2mag./rmag;
theta2 = acos(cost2);

% see which solutions satisfy Snell's law
test = n1*sint1 - n2*sint2;
ind = find(abs(test)<1e-11);

% ...of these, which has the smallest electrical length
eleclength = n1*v1mag(ind)  + n2*v2mag(ind);
[val ind2] = min(eleclength);
ind = ind(ind2);

% if a solution remains, check if theta1 is less than pi/2
if ~isempty(ind) && (theta1(ind) <= pi/2)
    rx = rx(ind);
    ry = ry(ind);
    theta1 = theta1(ind);
    theta2 = theta2(ind);
    thetal = atan2(rx,y1-ry);
    v1mag = v1mag(ind);
    v2mag = v2mag(ind);
    
    % rotate solution back
    rrot = R'*[rx; ry];
    rx = rrot(1);
    ry = rrot(2);
else
    rx = nan;
    ry = nan;
    theta1 = nan;
    theta2 = nan;
    thetal = nan;
    v1mag = nan;
    v2mag = nan;
end



