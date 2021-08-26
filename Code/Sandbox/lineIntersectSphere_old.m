function yn = lineIntersectSphere(x1,x2,r);
% returns 1 if the line intesects the sphere, 0 otherwise
%
% x1,x2  Nx3 arrays of coordinate points
% r      radius of sphere
%
% yn     1 if line itersects, 0 if not

    l = (x2 - x1);
    l = l./sqrt(sum(l.^2,2));
    b = dot(l,x1,2).^2 - (dot(x1,x1,2) - r^2);
    yn = (b >= 0);
end
