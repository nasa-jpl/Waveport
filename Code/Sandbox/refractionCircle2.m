function [rx ry theta1 theta2 v1mag v2mag] = refractionCircle2(x_1,y_1,x_2,y_2,n_1,n_2,r)
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
% v1mag, v2mag  Lengths of the vectors from each point to the refraction point

if sqrt(x_1^2 + y_1^2) < r
    disp('exterior point is interior')
    return
end
if sqrt(x_2^2 + y_2^2) > r
    disp('interior point is exterior')
    return
end
n_1 = real(n_1);
n_2 = real(n_2);

% compute c coefficients
c_1 = n_1^2*x_1^2 - n_1^2*y_1^2;
c_2 = 2*n_1^2*x_1*y_1;
c_3 = -n_1^2*r^2*x_1^2;
c_4 = 2*x_1;
c_5 = 2*y_1;
c_6 = -r^2 - x_1^2 - y_1^2;
c_7 = n_2^2*x_2^2 - n_2^2*y_2^2;
c_8 = 2*n_2^2*x_2*y_2;
c_9 = -n_2^2*r^2*x_2^2;
c_10 = 2*x_2;
c_11 = 2*y_2;
c_12 = -r^2 - x_2^2 - y_2^2;

% compute b coefficients
b_1 = c_1*c_10 - c_4*c_7 - c_2*c_11 + c_5*c_8;
b_2 = c_1*c_12 - c_6*c_7;
b_3 = c_3*c_10 - c_4*c_9 + (c_2*c_11 - c_5*c_8)*r^2;
b_4 = c_3*c_12 - c_6*c_9;
b_5 = c_4*c_8 + c_5*c_7 - c_1*c_11 - c_2*c_10;
b_6 = c_6*c_8 - c_2*c_12;
b_7 = c_5*c_9 - c_3*c_11;

% compute a coefficients
a_6 = b_1^2 + b_5^2;
a_5 = 2*(b_1*b_2 + b_5*b_6);
a_4 = b_2^2 - b_5^2*r^2 + 2*b_7*b_5 + b_6^2 + 2*b_1*b_3;
a_3 = 2*(-b_5*b_6*r^2 + b_1*b_4 + b_2*b_3 + b_6*b_7);
a_2 = 2*b_2*b_4 - r^2*(b_6^2 + 2*b_5*b_7) + b_3^2 + b_7^2;
a_1 = 2*(-b_6*b_7*r^2 + b_3*b_4);
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
r1 = repmat([x_1 y_1 0],Nr,1);
r2 = repmat([x_2 y_2 0],Nr,1);
v1 = r1 - rr;
v2 = r2 - rr;
rmag = sqrt(sum(rr.^2,2));
v1mag = sqrt(sum(v1.^2,2));
v2mag = sqrt(sum(v2.^2,2));
sint1 = cross(rr,v1,2)./v1mag./rmag;
sint2 = cross(-rr,v2,2)./v2mag./rmag;
sint1 = abs(sint1(:,3));
sint2 = abs(sint2(:,3));

% compute theta1 and theta2 via dot product
cost1 = dot(rr,v1,2)./v1mag./rmag;
theta1 = acos(cost1);
cost2 = dot(-rr,v2,2)./v2mag./rmag;
theta2 = acos(cost2);

% see which solutions satisfy Snell's law
test = n_1*sint1 - n_2*sint2;
test = n_1*sin(theta1) - n_2*sin(theta2);

ind = find(abs(test)<1e-11);

% ...of these, which has the smallest electrical length
eleclength = n_1*v1mag(ind)  + n_2*v2mag(ind);
[val ind2] = min(eleclength);
ind = ind(ind2);

% if a solution remains, check if theta1 is less than pi/2
if ~isempty(ind) && (theta1(ind) <= pi/2)
    rx = rx(ind);
    ry = ry(ind);
    theta1 = theta1(ind);
    theta2 = theta2(ind);
    v1mag = v1mag(ind);
    v2mag = v2mag(ind);
else
    rx = nan;
    ry = nan;
    theta1 = nan;
    theta2 = nan;
    v1mag = nan;
    v2mag = nan;
end



