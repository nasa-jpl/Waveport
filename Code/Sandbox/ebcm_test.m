



%% ebcm surface discretization

% define surface points on a sphere
a = 1;
[th phi] = discoball(30);
[x y z] = sph2cart(a*ones(size(th)),th,phi);

% compute Delaunay triangulation
dt = delaunayTriangulation(x,y,z);

% take the convex hull of these points
[ch v] = convexHull(dt);

% plot
if 0
    trisurf(ch,dt.Points(:,1),dt.Points(:,2),dt.Points(:,3))
end

% extract three points of each triangle
r1 = [x(ch(:,1)) y(ch(:,1)) z(ch(:,1))];
r2 = [x(ch(:,2)) y(ch(:,2)) z(ch(:,2))];
r3 = [x(ch(:,3)) y(ch(:,3)) z(ch(:,3))];

% compute centroid of triangles, these will be the integration points
cen = (1/3)*(r1 + r2 + r3);

X = cen(:,1);
Y = cen(:,2);
Z = cen(:,3);
N = length(X);

% compute the two edge vectors from point 1
r21 = r2-r1;
r31 = r3-r1;

% take cross product of edge vectors, needed for area and normal vector
r21crossr31 = cross(r21,r31,2);

% area of each triangle, one half the cross product magnitude
dS = (1/2)*sqrt(sum(r21crossr31.^2,2));

% outward surface normal unit vector
nhat = unit(r21crossr31);
nx = nhat(:,1);
ny = nhat(:,2);
nz = nhat(:,3);

% plot
figure(2),clf,hold all
trisurf(ch,dt.Points(:,1),dt.Points(:,2),dt.Points(:,3))
plot3(cen(:,1),cen(:,2),cen(:,3),'.')
quiver3(cen(:,1),cen(:,2),cen(:,3),nhat(:,1),nhat(:,2),nhat(:,3))
hold off
axis equal
view([1 1 1])
grid minor


%% compute T-matrix via ebcm

% wavenumber in region 1
lam = 1;
k1 = 2*pi/lam;

% wavenumber in region 2
er2 = 2;
k2 = 2*pi/lam*sqrt(er2);

% maximum degree harmonic
L = 12;

% ebcm
[Tmm Tmn Tnm Tnn] = ebcm(L,X,Y,Z,dS,nx,ny,nz,k1,k2);

% assemble the full T-matrix
Tebcm = [Tmm Tmn; Tnm Tnn];

% plot
figure(3),clf
imagesc(abs(Tebcm))


%% compare ebcm to dielectric sphere

[Tmm, Tnn] = tmatrixDielectricSphere(L,a,k1,k2,mu1,mu2);
Tsphere = [Tmm; Tnn];
T = diag(Tebcm);

plot([real([Tsphere T])])

plot(imag([Tsphere T]))

plot(abs([Tsphere T]))


%%

% 
% Q(indM,indM) = Q(indM,indM) + const*M1hat_cross_rgN2 + N1hat_cross_rgM2; 
% Q(indM,indN) = Q(indM,indN) + const*M1hat_cross_rgM2 + N1hat_cross_rgN2;
% Q(indN,indM) = Q(indN,indM) + const*N1hat_cross_rgN2 + M1hat_cross_rgM2;
% Q(indN,indN) = Q(indN,indN) + const*N1hat_cross_rgM2 + M1hat_cross_rgN2;
% 
% RgQ(indM,indM) = RgQ(indM,indM) + const*rgM1hat_cross_rgN2 + rgN1hat_cross_rgM2;
% RgQ(indM,indN) = RgQ(indM,indN) + const*rgM1hat_cross_rgM2 + rgN1hat_cross_rgN2;
% RgQ(indN,indM) = RgQ(indN,indM) + const*rgN1hat_cross_rgN2 + rgM1hat_cross_rgM2;
% RgQ(indN,indN) = RgQ(indN,indN) + const*rgN1hat_cross_rgM2 + rgM1hat_cross_rgN2;
% 


 