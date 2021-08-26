

%% Tmatrix of dieletric sphere

L = 30;
a = 3;
k1 = 2*pi;
k2 = 2*pi*sqrt(2);

[Tmm, Tnn] = tmatrixDieletricSphere(L,a,k1,k2);

Lmax = 1.5*k1*a;
Nmax = Lmax^2 + 2*Lmax;

plot(abs([Tmm Tnn]))
myplot('T-Matrix Dieletric Sphere','(l,m) Linear Index','abs(T)')
leg=legend('Tmm','Tnn')



%% Tmatrix of PEC sphere


L = 30;
a = 3;
k = 2*pi;

[Tmm, Tnn] = tmatrixPECSphere(L,a,k);

Lmax = 1.2*k*a;
Nmax = Lmax^2 + 2*Lmax;

plot(abs([Tmm Tnn]))
myplot('T-Matrix PEC Sphere','(l,m) Linear Index','abs(T)')
leg=legend('Tmm','Tnn')



%% ebcm

%%% ebcm surface discretization

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


%%% compute T-matrix via ebcm

% wavenumber in region 1
lam = 1;
k1 = 2*pi/lam;

% wavenumber in region 2
er2 = 2;
k2 = 2*pi/lam*sqrt(er2);

% maximum degree harmonic
L = 10;

% ebcm
[Tmm Tmn Tnm Tnn] = ebcm(L,X,Y,Z,dS,nx,ny,nz,k1,k2);

% assemble the full T-matrix
Tebcm = [Tmm Tmn; Tnm Tnn];

% take diagonal 
T1 = diag(Tmm);
T2 = diag(Tnn);

% plot
figure(3),clf
imagesc(abs(Tebcm))
myplot('T-matrix via EBCM','(l'',m'') linear index','(l,m) linear index')


%%% compare ebcm solution to a dielectric sphere

% Analytical solution for T-matrix of dielectric sphere 
[Tmm_ana, Tnn_ana] = tmatrixDielectricSphere(L,a,k1,k2);

% take diagonal from the ebcm solution

figure(3),clf,hold all
plot(real([Tmm_ana T1]))
plot(imag([Tmm_ana T1]))
hold off
myplot('T-matrix of Sphere - Tmm','(l,m) linear index','Real/Imag T-matrix Elements')
leg = legend('Real T-matrix - Analytical','Real T-matrix - EBCM',...
    'Imag T-matrix - Analytical','Imag T-matrix - EBCM','location','southeast');


figure(4),clf,hold all
plot(real([Tnn_ana T2]))
plot(imag([Tnn_ana T2]))
hold off
myplot('T-matrix of Sphere - Tnn','(l,m) linear index','Real/Imag T-matrix Elements')
leg = legend('Real T-matrix - Analytical','Real T-matrix - EBCM',...
    'Imag T-matrix - Analytical','Imag T-matrix - EBCM','location','southeast');




%% BRCS PEC sphere

k = 1;

a = logspace(-1,2,1000)';
Na = length(a);
brcsall = zeros(Na,1);
brcsallnorm = zeros(Na,1);

for n=1:Na,
    brcsall(n) = brcs_pec_sphere(k,a(n));
    brcsallnorm(n) = (1/(pi*a(n)^2))*brcs_pec_sphere(k,a(n));
end

figure(1),clf,
plot(10*log10(k*a),10*log10(brcsall))
myplot('RCS of a PEC Sphere','ka','\sigma_{PEC} (dB)')

h2 = figure(2);,clf,
loglog(k*a,brcsallnorm)
myplot('RCS of a PEC Sphere - Normalized','ka','\sigma_{PEC} /\pi a^2')
ylim([1e-3 10])




%% compute_T_from_S

% harmonics
L = 10;
k = 1;

% ficticious T-matrix
rng(1);
N = L^2 + 2*L;
sz = N;
Tmm = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);
Tmn = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);
Tnm = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);
Tnn = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);

% sample the S-matrix at the nodes of quadrature so we can compare to the FMM version
I = 2*L + 1;
J = L + 1;
phi = 2*pi*(0:(I-1))/I;
[muj, wj] = legpts(J);
theta = acos(muj);
[theta_s phi_s] = meshgrid(theta,phi);
[theta_i phi_i] = meshgrid(theta,phi);

% Compute S-matrix from T-matrix at the nodes of quadrature
[Stt, Stp, Spt, Spp] = compute_S_from_T(Tmm,Tmn,Tnm,Tnn,L,k,theta_s,phi_s,theta_i,phi_i);

% Compute S-matrix from T-matrix using the FMM version 
[Stt2, Stp2, Spt2, Spp2] = convert_T_to_S(Tmm,Tmn,Tnm,Tnn,k,I,J,L);

% compare
mean(abs(Stt(:) - Stt2(:)))
mean(abs(Stp(:) - Stp2(:)))
mean(abs(Spt(:) - Spt2(:)))
mean(abs(Spp(:) - Spp2(:)))


%% far-field T-matrix

% harmonics
L = 10;
k = 1;

% ficticious T-matrix
rng(1);
N = L^2 + 2*L;
sz = N;
Tmm = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);
Tmn = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);
Tnm = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);
Tnn = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);

% Compute Tfar from T
[Tbb, Tbc, Tcb, Tcc]        = convert_T_to_Tfar(Tmm,Tmn,Tnm,Tnn,L,k);

% Compute T from T far
[Tmm2, Tmn2, Tnm2, Tnn2]    = convert_Tfar_to_T(Tbb,Tbc,Tcb,Tcc,L,k);

% Compare the conversion
mean(abs(Tmm(:) - Tmm2(:)))
mean(abs(Tmn(:) - Tmn2(:)))
mean(abs(Tnm(:) - Tnm2(:)))
mean(abs(Tnn(:) - Tnn2(:)))



%% compute_scs_from_tmatrix

% Test this on dielectric sphere

% Parameters for the sphere T-matrix
L = 30;
a = 3;
k1 = 2*pi;
k2 = 1.2*k1;
u1 = 1;
u2 = 1;

% Compute the T-matrix
[Tmm, Tnn] = tmatrixDieletricSphere(L,a,k1,k2,u1,u2);
Tmm = diag(Tmm);
Tmn = zeros(size(Tmm));
Tnm = zeros(size(Tmm));
Tnn = diag(Tnn);

% Define polarization angle, beta
beta = linspace(0,2*pi,10);

% Define incident directions
theta_i = linspace(0,pi,10);
phi_i = linspace(0,2*pi,10);

% Compute the SCS from the T-matrix
scs_sphere = compute_scs_from_tmatrix(Tmm,Tmn,Tnm,Tnn,L,k1,theta_i,phi_i,beta);

% Compute the SCS directly
scs_sphere2 = scs_dieletric_sphere(k1,k2,a);

%% compute_avescs_from_tmatrix

% Test this on a dielectric sphere

% Parameters for the sphere T-matrix
L = 30;
a = 3;
k1 = 2*pi;
k2 = 1.2*k1;
u1 = 1;
u2 = 1;

% Compute the T-matrix
[Tmm, Tnn] = tmatrixDieletricSphere(L,a,k1,k2,u1,u2);
Tmm = diag(Tmm);
Tmn = zeros(size(Tmm));
Tnm = zeros(size(Tmm));
Tnn = diag(Tnn);

% Compute the polarization averaged SCS from T-matrix
avescs = compute_avescs_from_tmatrix(Tmm,Tmn,Tnm,Tnn,k1);

% Compute the SCS for the sphere directly 
scs_sphere2 = scs_dieletric_sphere(k1,k2,a);







