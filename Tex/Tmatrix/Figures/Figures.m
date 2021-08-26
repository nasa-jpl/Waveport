

direc = '/Users/mshaynes/Desktop/Book/Waveport/Tex/Tmatrix/Figures/';


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
h2=figure(2),clf,hold all
trisurf(ch,dt.Points(:,1),dt.Points(:,2),dt.Points(:,3))
plot3(cen(:,1),cen(:,2),cen(:,3),'.')
quiver3(cen(:,1),cen(:,2),cen(:,3),nhat(:,1),nhat(:,2),nhat(:,3))
hold off
axis equal
view([1 1 1])
grid minor
axis(1.1*[-1 1 -1 1 -1 1])
title('EBCM Surface Discretization of a Dielectric Sphere')
xlabel('x (\lambda)')
ylabel('y (\lambda)')
zlabel('z (\lambda)')
set(gcf,'color','w')
set(gca,'fontsize',14)

if 0
    saveimage(h2,[direc 'ebcmsphere'],'epsc');
end

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
[Tmm_ana, Tnn_ana] = tmatrixDielectricSphere(L,a,k1,k2,mu1,mu2);

% take diagonal from the ebcm solution

h3 = figure(3),clf,hold all
plot(real([Tmm_ana T1]))
plot(imag([Tmm_ana T1]))
hold off
myplot('T-matrix of a Dielectric Sphere - Tmm','(l,m) linear index','Real/Imag T-matrix Elements')
leg = legend('Real T-matrix - Analytical','Real T-matrix - EBCM',...
    'Imag T-matrix - Analytical','Imag T-matrix - EBCM','location','southeast');


h4 = figure(4),clf,hold all
plot(real([Tnn_ana T2]))
plot(imag([Tnn_ana T2]))
hold off
myplot('T-matrix of a Dielectric Sphere - Tnn','(l,m) linear index','Real/Imag T-matrix Elements')
leg = legend('Real T-matrix - Analytical','Real T-matrix - EBCM',...
    'Imag T-matrix - Analytical','Imag T-matrix - EBCM','location','southeast');


if 0
    saveimage(h3,[direc 'ebcmTmm'],'epsc');
    saveimage(h4,[direc 'ebcmTnn'],'epsc');
end





%% BRCS PEC

k = 1;

a = logspace(-1,2,1000)';
Na = length(a);
brcsall = zeros(Na,1);
brcsallnorm = zeros(Na,1);

for n=1:Na,
    brcsall(n) = brcs_pec_sphere(k,a(n));
    brcsallnorm(n) = (1/(pi*a(n)^2))*brcsall(n);
end

figure(1),clf,
plot(10*log10(k*a),10*log10(brcsall))
myplot('Backscatter RCS PEC Sphere','ka','\sigma_{PEC}')

h2 = figure(2);,clf,
loglog(k*a,brcsallnorm)
myplot('Backscatter RCS of a PEC Sphere','ka','\sigma_{PEC} /\pi a^2')
ylim([1e-3 10])

if 0
    set(h2,'PaperPositionMode','auto')
    saveas(h2,[direc 'brcspecsphere'],'epsc');
end


%% BRCS dieletric sphere


a = logspace(-1,1,1000)';
a = linspace(0.01,3,30000)';
Na = length(a);
brcsall = zeros(Na,1);
brcsallnorm = zeros(Na,1);

er1 = 1;
er2 = 2.5;
lam = 1;
k1 = 2*pi/lam*sqrt(er1);
k2 = 2*pi/lam*sqrt(er2);

for n=1:Na,
    brcsallnorm(n) = (1/(pi*a(n)^2))*brcs_dielectric_sphere(k1,k2,a(n));
end

h2 = figure(2);,clf,
loglog(k*a,brcsallnorm)
myplot('Backscatter RCS of a Dieletric Sphere','ka','\sigma /\pi a^2')
ylim([1e-3 10])

h2 = figure(2);,clf,
semilogy(a/lam,(brcsallnorm))
myplot('Backscatter RCS of a Dieletric Sphere','a/\lambda_o','\sigma /\pi a^2')
ylim([0.01 100])
xlim([0 a(end)])

if 0
    set(h2,'PaperPositionMode','auto')
    saveas(h2,[direc 'brcsdieletricsphere'],'epsc');
end


%% BRCS dieletric sphere grid

a = linspace(0.01,2,2000)';
Na = length(a);

er1 = 1;
lam = 1;
k1 = 2*pi/lam*sqrt(er1);

er2 = linspace(1,5,400);
Ner = length(er2);
brcsallnorm = zeros(Na,Ner);

for m=1:Ner,
    mydisp(m,Ner)
    k2 = 2*pi/lam*sqrt(er2(m));
    for n=1:Na,
        brcsallnorm(n,m) = (1/(pi*a(n)^2))*brcs_dielectric_sphere(k1,k2,a(n));
    end
end

h2 = figure(1);,clf,
imagesc(a/lam,er2,10*log10(brcsallnorm)')
set(gca,'ydir','normal')
mycbar({'\sigma/\pia^2';'(dB)'},[-30 15])
myplot('Backscatter RCS of a Dieletric Sphere','a/\lambda_o','\epsilon_r')


if 0
    set(h2,'PaperPositionMode','auto')
    saveas(h2,[direc 'brcsdieletricspheregrid'],'epsc');
end


%% SCS PEC

k = 1;
a = logspace(-1,2,1000)';
Na = length(a);
scssall = zeros(Na,1);
scsallnorm = zeros(Na,1);

for n=1:Na,
    scsall(n) = scs_pec_sphere(k,a(n));
    scsallnorm(n) = (1/(pi*a(n)^2))*scsall(n);
end

figure(1),clf,
plot(10*log10(k*a),10*log10(scsall))
myplot('Scattering Cross Section of a PEC Sphere','ka','\sigma_{PEC}')

h2 = figure(2);,clf,
loglog(k*a,scsallnorm)
myplot('Scattering Cross Section of a PEC Sphere','ka','\sigma_{PEC} /\pi a^2')
ylim([1e-3 10])
%xlim([0.5 15])

if 0
    set(h2,'PaperPositionMode','auto')
    saveas(h2,[direc 'scspecsphere'],'epsc');
end


%% SCS dieletric sphere grid

a = linspace(0.01,2,2000)';
Na = length(a);

er1 = 1;
lam = 1;
k1 = 2*pi/lam*sqrt(er1);

er2 = linspace(1,5,400);
Ner = length(er2);
scsallnorm = zeros(Na,Ner);

for m=1:Ner,
    mydisp(m,Ner)
    k2 = 2*pi/lam*sqrt(er2(m));
    for n=1:Na,
        scsallnorm(n,m) = (1/(pi*a(n)^2))*scs_dielectric_sphere(k1,k2,a(n));
    end
end

h2 = figure(1);,clf,
imagesc(a/lam,er2,10*log10(scsallnorm)')
set(gca,'ydir','normal')
myplot('Scattering Cross Section of a Dieletric Sphere','a/\lambda_o','\epsilon_r')
mycbar({'\sigma/\pia^2';'(dB)'},[-20 10])


if 0
    set(h2,'PaperPositionMode','auto')
    saveas(h2,[direc 'scsdieletricspheregrid'],'epsc');
end





