function [rx ry rz theta1 theta2 thetal v1mag v2mag] = refractionSphere(x1,y1,z1,x2,y2,z2,n1,n2,r,type,nit)
% Refraction through a spherical interface, assumes n1, n2 are real.
%
% x1,y1,z1      Coordinates of exterior point, r1, centered on circle
% x2,y2,z2      Coordinates of interior point, r2, centered on circle
% n1            Index of refraction for exterior
% n2            Index of refraction for interior
% r             Radius circle
% type          [optional] Defaults to direct solution. 
%               Use 'dir' for direct solution, 'it' for iterative. 
% nit           [optional] Number of iterations when using 'it' option 
%
% rx,ry,rz      Coordinates of the refraction point
% theta1        Incidence angle
% theta2        Transmission angle
% thetal        Look angle
% v1mag, v2mag  Lengths of the vectors from each point to the refraction point
%
% Dependencies: refractionCircle, refractionCircleIterative

sz = size(x1);
n1 = real(n1);
n2 = real(n2);

% columnize the inputs
r1 = [x1(:) y1(:) z1(:)];
r2 = [x2(:) y2(:) z2(:)];
N = length(x1(:));

% project the points in the r1-r2 plane
% v_hat
r1mag = sqrt(sum(r1.^2,2));
v_hat = r1./repmat(r1mag,1,3);
% w_hat
r2xvhat = cross(r2,v_hat,2);
% check for zero values (interoir and exterior point are radially aligned)
ind = find(sum(abs(r2xvhat),2)==0);
% in these cases, create an arbitrary cross product by permuting r1
tmp = [r1(ind,2) r1(ind,3) r1(ind,1)];
r2xvhat(ind,:) = tmp;
w_hat = r2xvhat./repmat(sqrt(sum(r2xvhat.^2,2)),1,3);
% u_hat
u_hat = cross(v_hat,w_hat,2);

% in the circle frame:
% y coordinte of exterior point
c1x = zeros(N,1);
c1y = r1mag;
% x coordinate of interior point
c2x = dot(r2,u_hat,2);
c2y = dot(r2,v_hat,2);

% initialize outputs
rx = zeros(N,1);
ry = zeros(N,1);
theta1 = zeros(N,1);
theta2 = zeros(N,1);
thetal = zeros(N,1);
v1mag = zeros(N,1);
v2mag = zeros(N,1);

if nargin >= 10 && strcmp(type,'dir') || nargin >= 10 && isempty(type) % loop over direct solution
    for n=1:N,
        [rxt ryt t1 t2 tl v1m v2m] = refractionCircle(c1x(n),c1y(n),c2x(n),c2y(n),n1,n2,r);
        rx(n) = rxt;
        ry(n) = ryt;
        theta1(n) = t1;
        theta2(n) = t2;
        thetal(n) = tl;
        v1mag(n) = v1m;
        v2mag(n) = v2m;
    end
elseif nargin >= 10 && strcmp(type,'it') % vectorized iterative solution
    if ~isempty(nit)
        [rx ry theta1 theta2 thetal v1mag v2mag] = refractionCircleIterative(c1x,c1y,c2x,c2y,n1,n2,r,nit);
    else
        [rx ry theta1 theta2 thetal v1mag v2mag] = refractionCircleIterative(c1x,c1y,c2x,c2y,n1,n2,r);
    end   
else
    error('bad inputs')
end

% project refraction point back to 3D
rr = repmat(rx,1,3).*u_hat + repmat(ry,1,3).*v_hat;

% reshape to original size of inputs
rx=reshape(rr(:,1),sz);
ry=reshape(rr(:,2),sz);
rz=reshape(rr(:,3),sz);
theta1=reshape(theta1,sz);
theta2=reshape(theta2,sz);
thetal=reshape(thetal,sz);
v1mag=reshape(v1mag,sz);
v2mag=reshape(v2mag,sz);