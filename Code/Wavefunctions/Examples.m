
%% Spherical harmonic indexing

% linear index to harmonics
ind = (1:8)';
% vector waves
[l m] = ind2lm(ind)
% scalar waves
[l m] = ind2lm(ind,'mono')

% harmonics to linear index
l = 4;
m = -3;
% vector waves
ind = lm2ind(l,m)
% scalar waves
ind = lm2ind(l,m,'mono')

% (l,m) table
L = 5;
% vector waves
[tab] = lmtable(L)
% scalar waves
[tab] = lmtable(L,'mono')



%% sphericalY

%%% compute Ylm (no monopole)
% theta/phi grid
Nth = 30;
theta = linspace(0,pi,Nth);
phi = linspace(0,2*pi,2*Nth);
[T P] = ndgrid(theta,phi);

% maximum harmonic, L
L = 3;
% total number of harmonics for L
N = L^2 + 2*L;

% compute spherical harmonics
ylm = sphericalY(L,T,P);


%%% plot individual harmonics
l = 3;
m = 2;
% get linear index for this harmonic 
ind = lm2ind(l,m);
% grab that harmonic and reshape into theta/phi grid
ylm_ind = reshape(ylm(:,ind),size(T));

figure(1),clf,
imagesc(phi,theta,real(ylm_ind));
myplot(['Spherical Harmonic, (l,m) = (' num2str(l) ', ' num2str(m) '), Real'],...
    '\phi (rad)','\theta (rad)')

figure(2),clf,
imagesc(phi,theta,imag(ylm_ind));
myplot(['Spherical Harmonic, (l,m) = (' num2str(l) ', ' num2str(m) '), Imag'],...
    '\phi (rad)','\theta (rad)')


%%% example of computing expansion sum
% create N = L^2 + 2*L random expansion coefficients through har
alm = rand(N,1)-0.5 + 1i*(rand(N,1)-0.5);

% compute expansion sum as a matrix vector multiply over coeffients
f = ylm*alm;

% reshape the result into the theta/phi grid
f = reshape(f,size(T));

% plot
figure(3),clf
imagesc(phi,theta,real(f));
myplot(['Spherical Function, Real'],'\phi (rad)','\theta (rad)')

figure(4),clf,
imagesc(phi,theta,imag(f));
myplot(['Spherical Function, Imag'],'\phi (rad)','\theta (rad)')


%%% test the spherical harmonic normalization by computing the numerical
%%% intergral over the sphere of the dot product
% create a fine theta/phi grid
Nth = 200;
theta = linspace(0,pi,Nth);
phi = linspace(0,2*pi,2*Nth);
[T P] = ndgrid(theta,phi);
L = 3;
N = L^2 + 2*L;

% compute all harmonics (no monopole)
ylm = sphericalY(L,T,P);

% compute step
dt = theta(2)-theta(1);
dp = phi(2)-phi(1);

% compute the dot product and integrate over the 
% sphere of each harmonic at once. Result is close to 1 and gets better
% with more points.
dt*dp*sum(ylm.*conj(ylm).*repmat(sin(T(:)),1,N),1)


%% psilm, scalar wave functions

% create r/theta/phi grid for two radial points
r = [2 4];
Nth = 40;
theta = linspace(0,pi,Nth);
phi = linspace(0,2*pi,2*Nth);
[R T P] = ndgrid(r,theta,phi);
L = 5;
N = L^2 + 2*L + 1;  % total number of scalar harmonics, including monopole 
k = 1; % background wave number equal to 1

% compute scalar wave functions (includes monopole)
psi = psilm(L,k,R,T,P);


%%% plot one harmonic for one radial point over the sphere
l = 4;
m = 4;
ind = lm2ind(l,m,'mono'); % linear index of the harmonic

% first grab the harmoic and reshape into the r/theta/phi grid
psi2 = reshape(psi(:,ind),size(T));

% second, grab the radial point
r_ind = 1; % look at first radial point
psi2 = squeeze(psi2(r_ind,:,:));

% plot
figure(1),clf,
imagesc(phi,theta,real(psi2)),colorbar


%%% example of computing expansion sum
% create N = L^2 + 2*L random expansion coefficients through har
alm = rand(N,1)-0.5 + 1i*(rand(N,1)-0.5);

% compute expansion sum as a matrix vector multiply over coeffients
f = psi*alm;

% reshape the result into the theta/phi grid
f = reshape(f,size(T));

% plot the second radial index
r_ind = 2;
f = squeeze(f(r_ind,:,:));

% plot
figure(3),clf
imagesc(phi,theta,real(f));
myplot(['Scalar Wave Function, Real'],...
    '\phi (rad)','\theta (rad)')

figure(4),clf,
imagesc(phi,theta,imag(f));
myplot(['Scalar Wave Function, Imag'],...
    '\phi (rad)','\theta (rad)')


%% scalarPlaneWaveCoef, test of scalar plane wave coefficients

% choose a propagation direction (\theta_k,\phi_k)
thetak = pi/4;
phik = pi/4;

% wave number equal to 1
k = 1; 

% convert to Cartersian unit wave vector direction
[kx ky kz] = sph2cart(1,thetak,phik);

% create a 3D grid of points around the origin
x = linspace(-3,3,20);
y = x;
z = x;
[X Y Z] = ndgrid(x,y,z);

% compute plane wave analytically over the grid
eikr = exp(1i*k*(kx*X + ky*Y + kz*Z));

% convert the grid points to spherical coordinates 
[R T P] = cart2sph(X,Y,Z);

% choose a maximum harmonic for the expansion
L = 30;

% compute the regular wave functions at each point on the grid
psi = psilm(L,k,R,T,P,'rg');

% compute the scalar plane wave coefficients
alm = scalarPlaneWaveCoef(L,thetak,phik);

% evaluate the expansion as a sum over wave functions and coefficients for
% each harmonic
f = psi*alm;

% plot the difference between the expansion and analytic.
% do real and imag at once, just vectorize all the points.
% change L and watch the error change.
% the points at the edges of the plots are furthest from the origin.
figure(2),clf,
hold all
plot(real(f(:) - eikr(:)))
plot(imag(f(:) - eikr(:)))
hold off
myplot({'Scalar Plane Wave Expansion Error';'Expansion - Analytical'},'Point','Error')
legend('Real','Imag')


%% BC, vector spherical harmonics

% create a theta/phi grid
Nth = 40;
theta = linspace(0,pi,Nth);
phi = linspace(0,2*pi,2*Nth);
[T P] = ndgrid(theta,phi);

% choose maximum harmonic degree
L = 5;

% total number of harmonics
N = L^2 + 2*L;

% compute the fully normallized vector spherical wave harmonics
[Bth Bphi Cth Cphi] = BC(L,T,P);

% plot one
l = 4;
m = 4;
ind = lm2ind(l,m);
type = 3;
switch type
    case 1
        tmp = reshape(Bth(:,ind),size(T));
    case 2
        tmp = reshape(Bphi(:,ind),size(T));
    case 3
        tmp = reshape(Cth(:,ind),size(T));
    case 4
        tmp = reshape(Cphi(:,ind),size(T));
end


figure(1),clf,
imagesc(phi,theta,real(tmp));
myplot(['Vector Spherical Harmonic, (l,m) = (' num2str(l) ', ' num2str(m) '), Real'],...
    '\phi (rad)','\theta (rad)')

figure(2),clf,
imagesc(phi,theta,imag(tmp));
myplot(['Vector Spherical Harmonic, (l,m) = (' num2str(l) ', ' num2str(m) '), Imag'],...
    '\phi (rad)','\theta (rad)')



%%% BCmult, test the multiplication helper function
% create random coefficients with N total harmonics 
blm = rand(N,1);
clm = rand(N,1);
% use the outputs from BC and compute the theta/phi vector field components
[Fth, Fphi] = BCmult(Bth,Bphi,Cth,Cphi,blm,clm);

% reshape and plot
type = 1;
switch type
    case 1
        F = reshape(Fth,size(T));
        strtype = 'F_{\theta}';
    case 2
        F = reshape(Fphi,size(T));
        strtype = 'F_{\phi}';
end

figure(1),clf
imagesc(phi,theta,real(F))
myplot([strtype ', Real'],'\phi (rad)','\theta (rad)')

figure(2),clf
imagesc(phi,theta,imag(F))
myplot([strtype ', Imag'],'\phi (rad)','\theta (rad)')


%%% Test orthonormality via numerical integration
L = 5;
N = L^2 + 2*L;
sampletype = 2;
switch sampletype
    case 1  % trapezoidal integration 
        Nth = 200;
        theta = linspace(0,pi,Nth);
        phi = linspace(0,2*pi,2*Nth);
        [T P] = ndgrid(theta,phi);
        dt = theta(2)-theta(1);
        dp = phi(2)-phi(1);
        const = dt*dp;
        W = sin(T); 
    case 2 % Gaussian quadrature (see other chapters about this)
        I = 2*L+1;
        J = L+1;
        [xj wj] = legpts(J);
        theta = acos(xj);
        phi = (0:(I-1))*2*pi/I;
        [T P] = ndgrid(theta,phi);
        [W ~] = ndgrid(wj,phi);
        const = 2*pi/(2*L+1);
end
[Bth Bphi Cth Cphi] = BC(L,T,P,'norm');

% Test the normalization for B and C for each self term harmonic
const*sum((Bth.*conj(Bth) + Bphi.*conj(Bphi)).*repmat(W(:),1,N),1)
const*sum((Cth.*conj(Cth) + Cphi.*conj(Cphi)).*repmat(W(:),1,N),1)
    
% Test orthonormality of Bth, Bphi across all harmonics
mat = zeros(N,N);
for n=1:N,
for m=1:N,
    mat(n,m) = const*sum((Bth(:,n).*conj(Bth(:,m)) + Bphi(:,n).*conj(Bphi(:,m))).*W(:),1);
end
end
figure(1),clf
imagesc(real(mat)),colorbar
myplot('Orthonormalization, Real Part','B_{l,m}(\theta,\phi)','B_{l,m}^*(\theta,\phi)')
figure(2),clf
imagesc(imag(mat)),colorbar
myplot('Orthonormalization, Imag Part','B_{l,m}(\theta,\phi)','B_{l,m}^*(\theta,\phi)')

% Test orthonormality of Cth, Cphi across all harmonics
mat = zeros(N,N);
for n=1:N,
for m=1:N,
    mat(n,m) = const*sum((Cth(:,n).*conj(Cth(:,m)) + Cphi(:,n).*conj(Cphi(:,m))).*W(:),1);
end
end
figure(1),clf
imagesc(real(mat)),colorbar
myplot('Orthonormalization, Real Part','C_{l,m}(\theta,\phi)','C_{l,m}^*(\theta,\phi)')
figure(2),clf
imagesc(imag(mat)),colorbar
myplot('Orthonormalization, Imag Part','C_{l,m}(\theta,\phi)','C_{l,m}^*(\theta,\phi)')

% Test orthogonatlity of B and C across all harmonics (works best with quadrature sampling)
mat = zeros(N,N);
for n=1:N,
for m=1:N,
    mat(n,m) = const*sum((Bth(:,n).*conj(Cth(:,m)) + Bphi(:,n).*conj(Cphi(:,m))).*W(:),1);
end
end
figure(1),clf
imagesc(real(mat)),colorbar
myplot('Orthonormalization, Real Part','B_{l,m}(\theta,\phi)','C_{l,m}^*(\theta,\phi)')
figure(2),clf
imagesc(imag(mat)),colorbar
myplot('Orthonormalization, Imag Part','B_{l,m}(\theta,\phi)','C_{l,m}^*(\theta,\phi)')




%% MN, vector spherical wave functions

Nth = 20;
r = [2 4];
theta = linspace(0,pi,Nth);
phi = linspace(0,2*pi,2*Nth);
[R T P] = ndgrid(r,theta,phi);
L = 5;
N = L^2 + 2*L;
k = 1;

% set the strings
rgstr = 'rg';
hatstr = 'hat';
normstr = 'norm';
% compute vector spherical wave functions
[Mth Mphi Nr Nth Nphi] = MN(L,k,R,T,P,rgstr,hatstr,normstr);

% create random expansion coefficients for N total harmonics
alm = rand(N,1);
blm = rand(N,1);
% use outputs from MN and coefficients to compute the r/theta/phi electric field components
[Er Eth Ephi] = MNmult(Mth,Mphi,Nr,Nth,Nphi,alm,blm);

% reshape one component
type = 3;
switch type
    case 1
        F = reshape(Er,size(T));
        strtype = 'E_{r}';
    case 2
        F = reshape(Eth,size(T));
        strtype = 'E_{\theta}';
    case 3
        F = reshape(Ephi,size(T));
        strtype = 'E_{\phi}';
end
% grab one radial position
r_ind = 1;
F = squeeze(F(r_ind,:,:));

% plot
figure(1),clf
imagesc(phi,theta,real(F)),colorbar
myplot([strtype ', Real'],'\phi (rad)','\theta (rad)')

figure(2),clf
imagesc(phi,theta,imag(F)),colorbar
myplot([strtype ', Imag'],'\phi (rad)','\theta (rad)')


%% Test MN cross product relations

% wavenumber and radius
k = 1;
r = 10;

%%% Test orthonormality via numerical integration
L = 5;
N = L^2 + 2*L;
sampletype = 2;
tab = lmtable(L);
l = tab(:,1);
switch sampletype
    case 1  % trapezoidal integration 
        Nth = 200;
        theta = linspace(0,pi,Nth);
        phi = linspace(0,2*pi,2*Nth);
        [T P] = ndgrid(theta,phi);
        dt = theta(2)-theta(1);
        dp = phi(2)-phi(1);
        const = dt*dp;
        W = sin(T); 
    case 2 % Gaussian quadrature (see other chapters about this)
        I = 2*L+1;
        J = L+1;
        [xj wj] = legpts(J);
        theta = acos(xj);
        phi = (0:(I-1))*2*pi/I;
        [T P] = ndgrid(theta,phi);
        [W ~] = ndgrid(wj,phi);
        const = 2*pi/(2*L+1);
end

% constant radius
R = r*ones(size(T));

% set the strings
rgstr = [];
hatstr = [];
normstr = 'norm';
% compute vector spherical wave functions
[Mth Mphi Nr Nth Nphi] = MN(L,k,R,T,P,rgstr,hatstr,normstr);
Mr = zeros(length(Mth(:,1)),1);

% set the strings
rgstr = [];
hatstr = 'hat';
normstr = 'norm';
% compute vector spherical wave functions
[Mth_hat Mphi_hat Nr_hat Nth_hat Nphi_hat] = MN(L,k,R,T,P,rgstr,hatstr,normstr);
Mr_hat = Mr;


% Test orthonormality of M and N_hat across all harmonics on sphere
mat = zeros(N,N);
for n=1:N,
for m=1:N,
    CC = cross([Mr Mth(:,n) Mphi(:,n)],[Nr_hat(:,m) Nth_hat(:,m) Nphi_hat(:,m)],2);
    r_dot_CC = CC(:,1);
    mat(n,m) = const*sum(r_dot_CC.*W(:),1);
end
end

% Test orthonormality of N cross M_hat across all harmonics on sphere
mat = zeros(N,N);
for n=1:N,
for m=1:N,
    CC = cross([Nr(:,n) Nth(:,n) Nphi(:,n)],[Mr_hat Mth_hat(:,m) Mphi_hat(:,m)],2);
    r_dot_CC = CC(:,1);
    mat(n,m) = const*sum(r_dot_CC.*W(:),1);
end
end

% Test orthonormality of M_hat cross N across all harmonics on sphere
mat = zeros(N,N);
for n=1:N,
for m=1:N,
    CC = cross([Mr_hat Mth_hat(:,n) Mphi_hat(:,n)],[Nr(:,m) Nth(:,m) Nphi(:,m)],2);
    r_dot_CC = CC(:,1);
    mat(n,m) = const*sum(r_dot_CC.*W(:),1);
end
end

% Test orthonormality of N_hat cross M across all harmonics on sphere
mat = zeros(N,N);
for n=1:N,
for m=1:N,
    CC = cross([Nr_hat(:,n) Nth_hat(:,n) Nphi_hat(:,n)],[Mr Mth(:,m) Mphi(:,m)],2);
    r_dot_CC = CC(:,1);
    mat(n,m) = const*sum(r_dot_CC.*W(:),1);
end
end


figure(1),clf
imagesc(real(mat)),colorbar
myplot('Orthonormalization, Real Part','M_{l,m}(\theta,\phi)','N_{l,m}^*(\theta,\phi)')
figure(2),clf
imagesc(imag(mat)),colorbar
myplot('Orthonormalization, Imag Part','M_{l,m}(\theta,\phi)','N_{l,m}^*(\theta,\phi)')


% check values of the cross product
% grab diagonal 
D = diag(mat);
D = D(:);

% compare to the analytical (this combination needs to go with M cross N_hat)
hh = sbesselh(l,k*r);
hp = sbesselhp2(l,k*r)/(k*r);
ana = hh(:).*hp(:);

figure(1),clf
plot([real([ana-D]) imag([ana-D])])
myplot('Cross Product Orthonormalization','(l,m) diagonal linear index)','Real/Imag Numerial and Analytic')


% Test cross product sign of N cross M for equal harmonics
mat = zeros(N,1);
for m=1:N,
    N_cross_M = cross([Nr(:,m) Nth(:,m) Nphi(:,m)],[Mr Mth(:,m) Mphi(:,m)],2);
    M_cross_N = cross([Mr Mth(:,m) Mphi(:,m)],[Nr(:,m) Nth(:,m) Nphi(:,m)],2);
    mat(m) = sum(mean(abs(N_cross_M(:) + M_cross_N(:))));
end

% Test cross product sign of N cross Mhat and Mhat cross N for equal harmonics
mat = zeros(N,1);
for m=1:N,
    N_cross_Mhat = cross([Nr(:,m) Nth(:,m) Nphi(:,m)],[Mr_hat Mth_hat(:,m) Mphi_hat(:,m)],2);
    Mhat_cross_N = cross([Mr_hat Mth_hat(:,m) Mphi_hat(:,m)],[Nr(:,m) Nth(:,m) Nphi(:,m)],2);
    mat(m) = sum(mean(abs(N_cross_Mhat(:) + Mhat_cross_N(:))));
end

% Test cross product sign of M cross Nhat and Nhat cross M for equal harmonics
mat = zeros(N,1);
for m=1:N,
    M_cross_Nhat = cross([Mr Mth(:,m) Mphi(:,m)],[Nr_hat(:,m) Nth_hat(:,m) Nphi_hat(:,m)],2);
    Nhat_cross_M = cross([Nr_hat(:,m) Nth_hat(:,m) Nphi_hat(:,m)],[Mr Mth(:,m) Mphi(:,m)],2);
    mat(m) = sum(mean(abs(M_cross_Nhat(:) + Nhat_cross_M(:))));
end

% Test cross product sign of M cross Nhat and N cross M_hat for equal harmonics
mat = zeros(N,1);
for m=1:N,
    M_cross_Nhat = cross([Mr Mth(:,m) Mphi(:,m)],[Nr_hat(:,m) Nth_hat(:,m) Nphi_hat(:,m)],2);
    N_cross_Mhat = cross([Nr(:,m) Nth(:,m) Nphi(:,m)],[Mr_hat Mth_hat(:,m) Mphi_hat(:,m)],2);
    mat(m) = sum(mean(abs(M_cross_Nhat(:) + N_cross_Mhat(:))));
end


figure(1),clf
imagesc(real(mat)),colorbar
myplot('Orthonormalization, Real Part','M_{l,m}(\theta,\phi)','N_{l,m}^*(\theta,\phi)')
figure(2),clf
imagesc(imag(mat)),colorbar
myplot('Orthonormalization, Imag Part','M_{l,m}(\theta,\phi)','N_{l,m}^*(\theta,\phi)')


% Test cross product as before across all harmonics
mat = zeros(N,1);
for n=1:N,
for m=1:N,
    N_cross_M = cross([Nr(:,n) Nth(:,n) Nphi(:,n)],[Mr Mth(:,m) Mphi(:,m)],2);
    M_cross_N = cross([Mr Mth(:,n) Mphi(:,n)],[Nr(:,m) Nth(:,m) Nphi(:,m)],2);
    mat(n,m) = sum(mean(abs(N_cross_M(:) + M_cross_N(:))));
end
end

mat = zeros(N,1);
for n=1:N,
for m=1:N,
    N_cross_Mhat = cross([Nr(:,n) Nth(:,n) Nphi(:,n)],[Mr_hat Mth_hat(:,m) Mphi_hat(:,m)],2);
    Mhat_cross_N = cross([Mr_hat Mth_hat(:,n) Mphi_hat(:,n)],[Nr(:,m) Nth(:,m) Nphi(:,m)],2);
    mat(n,m) = sum(mean(abs(N_cross_Mhat(:) + Mhat_cross_N(:))));
end
end

mat = zeros(N,1);
for n=1:N,
for m=1:N,
    M_cross_Nhat = cross([Mr Mth(:,n) Mphi(:,n)],[Nr_hat(:,m) Nth_hat(:,m) Nphi_hat(:,m)],2);
    Nhat_cross_M = cross([Nr_hat(:,n) Nth_hat(:,n) Nphi_hat(:,n)],[Mr Mth(:,m) Mphi(:,m)],2);
    mat(n,m) = sum(mean(abs(M_cross_Nhat(:) + Nhat_cross_M(:))));
end
end

figure(1),clf
imagesc(real(mat)),colorbar
myplot('Orthonormalization, Real Part','M_{l,m}(\theta,\phi)','N_{l,m}^*(\theta,\phi)')
figure(2),clf
imagesc(imag(mat)),colorbar
myplot('Orthonormalization, Imag Part','M_{l,m}(\theta,\phi)','N_{l,m}^*(\theta,\phi)')





%% vectorPlaneWaveCoef, vector plane wave coefficients 


% choose a propagation direction (\theta_k,\phi_k)
thetak = pi/3;
phik = pi/4;

% wave number equal to 1
k = 1; 

% convert to Cartersian unit wave vector direction
[kx ky kz] = sph2cart(1,thetak,phik);

% create E-field components as an arbitrary perpenduclar vector to k
z_hat = [0 0 1];
E = cross([kx ky kz],z_hat,2)
E = E/sqrt(sum(E.^2));

% create a 3D grid of points around the origin
x = linspace(-3,3,20);
y = x;
z = x;
[X Y Z] = ndgrid(x,y,z);

% compute the scalar plane wave analytically over the grid
eikr = exp(1i*k*(kx*X + ky*Y + kz*Z));

% compute analytic vector plane wave
Ex_ana = E(1)*eikr;
Ey_ana = E(2)*eikr;
Ez_ana = E(3)*eikr;

% convert the grid points to spherical coordinates 
[R T P] = cart2sph(X,Y,Z);

% choose a maximum harmonic for the expansion
L = 30;

% compute the regular vector wave functions at each point on the grid
% set the strings
rgstr = 'rg';
hatstr = [];
normstr = 'norm';

% compute vector spherical wave functions
[Mth Mphi Nr Nth Nphi] = MN(L,k,R,T,P,rgstr,hatstr,normstr);

% compute the vector plane wave coefficients
[alm blm] = vectorPlaneWaveCoef(L,E(1),E(2),E(3),thetak,phik);

% evaluate the expansion as a sum over wave functions and coefficients for
% each harmonic
[Er Eth Ephi] = MNmult(Mth,Mphi,Nr,Nth,Nphi,alm,blm);

% convert to Cartesian vector field
[Ex Ey Ez] = sph2cart(R(:),T(:),P(:),Er,Eth,Ephi);


% plot the difference between the expansion and analytic.
% do real and imag at once, just vectorize all the points.
% change L and watch the error change.
% the points at the edges of the plots are furthest from the origin.

% grab one component
type = 1;
switch type
    case 1
        F = Ex;
        F_ana = Ex_ana(:);
        strtype = 'E_{x}';
    case 2
        F = Ey;
        F_ana = Ey_ana(:);
        strtype = 'E_{y}';
    case 3
        F = Ez;
        F_ana = Ez_ana(:);
        strtype = 'E_{z}';
end

% plot the difference between the expansion and analytical
figure(2),clf,
hold all
plot(real(F - F_ana))
plot(imag(F - F_ana))
hold off
myplot({['Vector Plane Wave Expansion Error, '  strtype];'Expansion - Analytical'},'Point','Error')
legend('Real','Imag')


%% vectorPlaneWaveCoef, vector plane wave coefficients 


% choose a propagation direction (\theta_k,\phi_k)
thetak = pi/3;
phik = pi/4;

% wave number equal to 1
k = 1; 

% convert to Cartersian unit wave vector direction
[kx ky kz] = sph2cart(1,thetak,phik);

% create E-field components as an arbitrary perpenduclar vector to k
z_hat = [0 0 1];
E = cross([kx ky kz],z_hat,2)
E = E/sqrt(sum(E.^2));

% create a 3D grid of points around the origin
x = linspace(-3,3,20);
y = x;
z = x;
[X Y Z] = ndgrid(x,y,z);

% compute the scalar plane wave analytically over the grid
eikr = exp(1i*k*(kx*X + ky*Y + kz*Z));

% compute analytic vector plane wave
Ex_ana = E(1)*eikr;
Ey_ana = E(2)*eikr;
Ez_ana = E(3)*eikr;

% convert the grid points to spherical coordinates 
[R T P] = cart2sph(X,Y,Z);

% choose a maximum harmonic for the expansion
L = 30;

% compute the regular vector wave functions at each point on the grid
% set the strings
rgstr = 'rg';
hatstr = [];
normstr = 'norm';
% compute vector spherical wave functions
[Mth Mphi Nr Nth Nphi] = MN(L,k,R,T,P,rgstr,hatstr,normstr);

% compute the vector plane wave coefficients
[alm blm] = vectorPlaneWaveCoef(L,E(1),E(2),E(3),thetak,phik);

% evaluate the expansion as a sum over wave functions and coefficients for
% each harmonic
[Er Eth Ephi] = MNmult(Mth,Mphi,Nr,Nth,Nphi,alm,blm);

% convert to Cartesian vector field
[Ex Ey Ez] = sph2cart(R(:),T(:),P(:),Er,Eth,Ephi);


% plot the difference between the expansion and analytic.
% do real and imag at once, just vectorize all the points.
% change L and watch the error change.
% the points at the edges of the plots are furthest from the origin.

% grab one component
type = 1;
switch type
    case 1
        F = Ex;
        F_ana = Ex_ana(:);
        strtype = 'E_{x}';
    case 2
        F = Ey;
        F_ana = Ey_ana(:);
        strtype = 'E_{y}';
    case 3
        F = Ez;
        F_ana = Ez_ana(:);
        strtype = 'E_{z}';
end

% plot the difference between the expansion and analytical
figure(2),clf,
hold all
plot(real(F - F_ana))
plot(imag(F - F_ana))
hold off
myplot({['Vector Plane Wave Expansion Error, '  strtype];'Expansion - Analytical'},'Point','Error')
legend('Real','Imag')


%% vectorPlaneWaveCoefZ, vector plane wave coefficients for z-propagating plane wave

% wave number equal to 1
k = 1; 

% create E-field components for the plane wave (no z component)
E = [1 1 0];

% create a 3D grid of points around the origin
x = linspace(-3,3,20);
y = x;
z = x;
[X Y Z] = ndgrid(x,y,z);

% compute the scalar plane wave analytically over the grid in z-direction
eikr = exp(1i*k*Z);

% compute analytic vector plane wave
Ex_ana = E(1)*eikr;
Ey_ana = E(2)*eikr;
Ez_ana = E(3)*eikr;

% convert the grid points to spherical coordinates 
[R T P] = cart2sph(X,Y,Z);

% choose a maximum harmonic for the expansion
L = 30;

% compute the regular vector wave functions at each point on the grid
% set the strings
rgstr = 'rg';
hatstr = [];
normstr = 'norm';
% compute vector spherical wave functions
[Mth Mphi Nr Nth Nphi] = MN(L,k,R,T,P,rgstr,hatstr,normstr);

% compute the vector plane wave coefficients for z-propagating plane wave
[alm blm] = vectorPlaneWaveCoefZ(L,E(1),E(2));

% evaluate the expansion as a sum over wave functions and coefficients for
% each harmonic
[Er Eth Ephi] = MNmult(Mth,Mphi,Nr,Nth,Nphi,alm,blm);

% convert to Cartesian vector field
[Ex Ey Ez] = sph2cart(R(:),T(:),P(:),Er,Eth,Ephi);


% plot the difference between the expansion and analytic.
% do real and imag at once, just vectorize all the points.
% change L and watch the error change.
% the points at the edges of the plots are furthest from the origin.

% grab one component
type = 3;
switch type
    case 1
        F = Ex;
        F_ana = Ex_ana(:);
        strtype = 'E_{x}';
    case 2
        F = Ey;
        F_ana = Ey_ana(:);
        strtype = 'E_{y}';
    case 3
        F = Ez;
        F_ana = Ez_ana(:);
        strtype = 'E_{z}';
end

% plot the difference between the expansion and analytical
figure(2),clf,
hold all
plot(real(F - F_ana))
plot(imag(F - F_ana))
hold off
myplot({['Vector Plane Wave Expansion Error, '  strtype];'Expansion - Analytical'},'Point','Error')
legend('Real','Imag')

%%% can also check this with the previous routine
thetak = 0;
phik = 0;
[alm1 blm1] = vectorPlaneWaveCoef(L,E(1),E(2),E(3),0,0);
[alm2 blm2] = vectorPlaneWaveCoefZ(L,E(1),E(2));

% compute the means
mean(abs([alm1 - alm2]))
mean(abs([blm1 - blm2]))



