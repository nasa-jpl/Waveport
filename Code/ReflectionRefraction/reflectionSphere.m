function [rx ry rz] = reflectionSphere(x1,y1,z1,x2,y2,z2,r)
% Reflection point on a sphere centered at the origin
%
% x1,y1,z1  Cartesian positions of first point (any size)
% x2,y2,z2  Cartesian positions of second point (same size as x1,y1,z1)
% r:        Sphere radius
%
% rx,ry,rz: Vector to the reflection point on sphere
%
% Based on David Eberly, Geometric Tools, 2008
%
% Dependencies: lineIntersectSphere,quarticroots

% check that points are outside sphere
sz = size(x1);
r1 = [x1(:) y1(:) z1(:)];
r2 = [x2(:) y2(:) z2(:)];
N = length(x1(:));

% points that are internal to the sphere
isbad = or((sqrt(sum(r1.^2,2)) <= r),(sqrt(sum(r2.^2,2)) <= r));

% points for which the line intersects the sphere (therefore shadowed)
[~, ~, ~, ~, delta_ind] = lineIntersectSphere(r1,r2,r);
isintersect = zeros(N,1);
isintersect(delta_ind == 1) = 1;

% check to see if those points are in front or behind
isshaddow = and((sign(dot(r1,r2,2)) == -1),isintersect);

% normalize points to radius of sphere
L = r1/r;
S = r2/r;
a = dot(S,S,2);
b = dot(S,L,2);
c = dot(L,L,2);
% fill coefficient array
qa = 4*c.*(a.*c-b.^2);
qb = -4*(a.*c-b.^2);
qc = a+2*b+c-4*a.*c;
qd = 2*(a-b);
qe = a-1;
[x1 x2 x3 x4] = quarticroots(qa,qb,qc,qd,qe);
ybar = [x1 x2 x3 x4];

% solve for projection scale factors
xbar = (-2*c.*ybar.^2 + ybar + 1)./(2*b.*ybar+1);
u = zeros(N,2);
for n=1:4,
    ind = and(and(imag(ybar(:,n))==0,ybar(:,n)>0),xbar(:,n)>0);
    u(ind,1) = xbar(ind,n);
    u(ind,2) = ybar(ind,n);
end
xbar = u(:,1);
ybar = u(:,2);
        
% reproject
rr = r*(xbar.*S + ybar.*L);
rr(or(isbad,isshaddow),:) = nan;
rx = reshape(rr(:,1),sz);
ry = reshape(rr(:,2),sz);
rz = reshape(rr(:,3),sz);

