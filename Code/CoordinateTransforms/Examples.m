
% Test coordinate transformations

% Create a 3D array of Cartesian points and vector field
x = linspace(-1,1,10);
y = x;
z = x;
[X Y Z] = meshgrid(x,y,z);
Ax = -Y;
Ay = -Z;
Az = -X;


%%% Cartesian to Cylindrical points
% Cartesian points to cylindrical points
[rho phi z] = cart2cyl(X,Y,Z);

% Cylindrical points to Cartesian points
[x2 y2 z2] = cyl2cart(rho,phi,z);

% Check the average error 
mean(abs(x2(:) - X(:)))
mean(abs(y2(:) - Y(:)))
mean(abs(z2(:) - Z(:)))


%%% Cartesian to Cylindrical vector fields
% Cartesian vector field to cylindrical vector field
[Arho Aphi Az] = cart2cyl(X,Y,Z,Ax,Ay,Az);

% Cylindrical vector field to Cartesian vector field
[Ax2 Ay2 Az2] = cyl2cart(rho,phi,z,Arho,Aphi,Az);

% Check the average error 
mean(abs(Ax2(:) - Ax(:)))
mean(abs(Ay2(:) - Ay(:)))
mean(abs(Az2(:) - Az(:)))


%%% Cartesian to Spherical points
% Cartesian points to spherical points
[r th phi] = cart2sph(X,Y,Z);

% Spherical points to Cartesian points 
[x2 y2 z2] = sph2cart(r,th,phi);

% Check the average error 
mean(abs(x2(:) - X(:)))
mean(abs(y2(:) - Y(:)))
mean(abs(z2(:) - Z(:)))


%%% Cartesian to Spherical vector fields
% Cartesian vector field to spherical vector field
[Ar Ath Aphi] = cart2sph(X,Y,Z,Ax,Ay,Az);

% Spherical vector field to Cartesian vector field
[Ax2 Ay2 Az2] = sph2cart(r,th,phi,Ar,Ath,Aphi);

% Check the average error 
mean(abs(Ax2(:) - Ax(:)))
mean(abs(Ay2(:) - Ay(:)))
mean(abs(Az2(:) - Az(:)))


%%% Cylindrical to Spherical points
% Cylindrical to spherical points
[r th phi] = cyl2sph(X,Y,Z);

% Spherical to cylindrical points
[x2 y2 z2] = sph2cyl(r,th,phi);

% Check the average error 
mean(abs(x2(:) - X(:)))
mean(abs(y2(:) - Y(:)))
mean(abs(z2(:) - Z(:)))


%%% Cylindrical to Spherical vector fields
% Cylindrical to spherical vector fields
[Ar Ath Aphi] = cyl2sph(X,Y,Z,Ax,Ay,Az);

% Spherical to cylindrical vector fields
[Ax2 Ay2 Az2] = sph2cyl(r,th,phi,Ar,Ath,Aphi);

% Check the average error 
mean(abs(Ax2(:) - Ax(:)))
mean(abs(Ay2(:) - Ay(:)))
mean(abs(Az2(:) - Az(:)))


%% Helper functions

%%% version 1 
N = 10;
Ndim = 4;
v = rand(N,Ndim);

% test mag
m = mag(v);
sqrt(sum(v.^2,2))

% test unit
vhat = unit(v);
sqrt(sum(vhat.^2,2))

%%% version 2
vx = rand(N);
vy = rand(N);
vz = rand(N);

% test mag
m = mag(vx,vy,vz);
mt = sqrt(vx.^2 + vy.^2 + vz.^2);
mean(abs(m(:)-mt(:)))

% test unit
[vhatx vhaty vhatz] = unit(vx,vy,vz);
sqrt(vhatx.^2+vhaty.^2+vhatz.^2)






