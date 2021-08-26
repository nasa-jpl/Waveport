function [i1 i2 l1 l2 delta_ind] = lineIntersectSphere(r1,r2,a,c);
% Line-sphere intersection
%
% r1            1x3 or Nx3 tail of points describing the line
% r2            1x3 or Nx3 head of points paired with x1
% a             radius of sphere
% c             [optional] 1x3 center of sphere, default to the origin
%
% i1,i2         Nx3 intersection points
% l1,l2         Nx3 distances from x1 to two intersection points
% delta_ind     Nx1 Determinant indicator  
%               -1: no intersection
%               0:  line is tangent
%               1:  2 intersections

N1 = length(r1(:,1));
N2 = length(r2(:,1));
if N1 == 1
    r1 = repmat(r1,N2,1);
elseif N2 == 1
    r2 = repmat(r2,N1,1);
elseif N1 ~= N2
    error('bad inputs')
end
N = max([N1 N2]);

% if no center defined, default to orgin
if nargin < 4
    c = [0 0 0];
end
c = repmat(c,N,1);

% unit vector from x1 to x2
uhat  = (r2-r1)./sqrt(dot(r2-r1,r2-r1,2));

% determinant 
delta = dot(uhat,(r1-c),2).^2 - (dot(r1-c,r1-c,2)-a^2);
 
% indicators
delta_ind = zeros(N,1);
delta_ind(delta < 0) = -1;
delta_ind(delta == 0) = 0;
delta_ind(delta > 0) = 1;

% two solutions of line lengths
part = -dot(uhat,(r1-c),2);
l1 = part + sqrt(delta);
l2 = part - sqrt(delta);

% two solutions for intersection points
i1 = r1 + l1.*uhat;
i2 = r1 + l2.*uhat;



