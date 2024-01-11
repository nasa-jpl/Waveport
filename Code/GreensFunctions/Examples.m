


%% 2D far field Green's function, relative to zero

% wavelength = 1
k = 2*pi;
% define radial variable
rho  = linspace(0,10,500)';

% 2D Green's function exact
h = (i/4)*besselh(0,k*rho);

% Far-field 2D Green's function
hfar = sqrt(1./(8*pi*k*rho)).*exp(i*(k*rho+pi/4));

% Compare
figure(1),clf,hold all
plot(rho,[real(h) imag(h)])
plot(rho,[real(hfar) imag(hfar)])
hold off
legend('Real g(\rho,\rho)','Imag g(\rho,\rho',...
    'Real far g(\rho,\rho)','Imag far g(\rho,\rho')


%% 2D far field Green's function, relative to a point

% Define 2D grid
k = 2*pi;
x = linspace(-50,50,500);
y = linspace(-50,50,500);
[X Y] = meshgrid(x,y);

% point source location
xp = [0 3];

% radial variable
rho = sqrt((X-xp(1)).^2 + (Y-xp(2)).^2);

% 2D Green's function exact
g1 = (1i/4)*besselh(0,k*rho);

% Far-field 2D Green's function
rho2 = sqrt(X.^2 + Y.^2);
X_hat = X./rho2;
Y_hat = Y./rho2;
rhodot = xp(1)*X_hat+xp(2)*Y_hat;
g2 = sqrt(1./(8*pi*k*rho2)).*exp(1i*pi/4).*exp(1i*k*rho2).*exp(-1i*k*rhodot);

% Plot exact Green's function
figure(1),clf
imagesc(x,y,real(g1)),colorbar
axis equal
caxis([-0.05 0.05])
myplot('2D Green''s Function','x (\lambda)','y (\lambda)')

% Plot far-field result
figure(2),clf
imagesc(x,y,real(g2)),colorbar
axis equal
caxis([-0.05 0.05])
myplot('Far-field 2D Green''s Function','x (\lambda)','y (\lambda)')

% Plot difference
figure(3),clf
imagesc(x,y,real(g2-g1)),colorbar
axis equal
caxis([-0.05 0.05])
myplot('Difference','x (\lambda)','y (\lambda)')

%% 2D far field Green's function, relative to a point, far cut

% lambda = 1
k = 2*pi;

% define 1D cut
X = linspace(100,110,1000)';
Y = X;

% source location
xp = [0 3];

% radial variable
rho = sqrt((X-xp(1)).^2 + (Y-xp(2)).^2);

% 2D Green's function exact
g1 = (1i/4)*besselh(0,k*rho);

% Far-field 2D Green's function
rho2 = sqrt(X.^2 + Y.^2);
X_hat = X./rho2;
Y_hat = Y./rho2;
rhodot = xp(1)*X_hat+xp(2)*Y_hat;
g2 = sqrt(1./(8*pi*k*rho2)).*exp(1i*pi/4).*exp(1i*k*rho2).*exp(-1i*k*rhodot);

% Plot both
figure(1),clf
plot(X,[real([g1 g2]) imag([g1 g2])])
myplot('Exact and Far-field 2D Green''s Function','\rho (\lambda)','Read/Imag')

%% dyadicGreens

% wavelength = 1
k = 2*pi;

% define 2D cut in XY plane
x = linspace(-5,5,300);
y = linspace(-5,5,300);
[X Y] = meshgrid(x,y);
Z = zeros(size(X));

% compute unique components of the dyadic Green's function
[Gxx Gyy Gzz Gxy Gxz Gyz] = dyadicGreens(k,X,Y,Z);

% plot 
cx = 0.1*[-1 1];
figure(1),clf
subplot(321),imagesc(x,y,real(Gxx)),caxis(cx),axis equal
myplot('Gxx','x (\lambda)','y (\lambda)')

subplot(323),imagesc(x,y,real(Gyy)),caxis(cx),axis equal
myplot('Gyy','x (\lambda)','y (\lambda)')

subplot(325),imagesc(x,y,real(Gzz)),caxis(cx),axis equal
myplot('Gzz','x (\lambda)','y (\lambda)')

subplot(322),imagesc(x,y,real(Gxy)),caxis(cx),axis equal
myplot('Gxy','x (\lambda)','y (\lambda)')

subplot(324),imagesc(x,y,real(Gxz)),caxis(cx),axis equal
myplot('Gxz','x (\lambda)','y (\lambda)')

subplot(326),imagesc(x,y,real(Gyz)),caxis(cx),axis equal
myplot('Gyz','x (\lambda)','y (\lambda)')



%% curlDyadicGreens

% wavelength = 1
k = 2*pi;

% define 2D cut in XY plane
x = linspace(-5,5,300);
y = linspace(-5,5,300);
[X Y] = meshgrid(x,y);
Z = zeros(size(X));

% compute unique components of the dyadic Green's function
[cGxy cGyx cGxz cGzx cGyz cGzy] = curlDyadicGreens(k,X,Y,Z);

% plot 
cx = 0.5*[-1 1];
figure(1),clf
subplot(321),imagesc(x,y,real(cGxy)),caxis(cx),axis equal
myplot('[Curl G]xy','x (\lambda)','y (\lambda)')

subplot(323),imagesc(x,y,real(cGyx)),caxis(cx),axis equal
myplot('[Curl G]yx','x (\lambda)','y (\lambda)')

subplot(325),imagesc(x,y,real(cGxz)),caxis(cx),axis equal
myplot('[Curl G]xz','x (\lambda)','y (\lambda)')

subplot(322),imagesc(x,y,real(cGzx)),caxis(cx),axis equal
myplot('[Curl G]zx','x (\lambda)','y (\lambda)')

subplot(324),imagesc(x,y,real(cGyz)),caxis(cx),axis equal
myplot('[Curl G]yz','x (\lambda)','y (\lambda)')

subplot(326),imagesc(x,y,real(cGzy)),caxis(cx),axis equal
myplot('[Curl G]zy','x (\lambda)','y (\lambda)')


%% volintGreens

kb = 2*pi;
a = 1/10;
% str = '2D';
% str = '3D';
str = 'dyadic'
[sing delta] = volintGreens(a,kb,str)


%% momGmatrix2D

k = 2*pi;

% create a dense square grid of points
N = 10;
x = linspace(-2,2,N);
y = linspace(-2,2,N);
[X Y] = meshgrid(x,y);
dx = x(2)-x(1);

% return full matrix, set primed and unprimed coordinates to be the same
Xp = X;
Yp = Y;

% build the matrix
G = momGmatrix2D(dx,k,X,Y,Xp,Yp);

% plot
figure(1),clf
imagesc(real(G))
myplot('MoM Green''s Function Matrix - 2D','Point index (x'',y'')','Point index (x,y)')
colorbar
mycbar('real(G)',1e-2*[-1 1])

figure(2),clf
imagesc(imag(G))
myplot('MoM Green''s Function Matrix - 2D','Point index (x'',y'')','Point index (x,y)')
colorbar
mycbar('imag(G)',1e-2*[-1 1])

figure(3),clf
imagesc(abs(G))
myplot('MoM Green''s Function Matrix - 2D','Point index (x'',y'')','Point index (x,y)')
colorbar
mycbar('abs(G)',[0 1e-2])


%% momGmatrix3D

k = 2*pi;

% create a dense square grid of points
N = 10;
x = linspace(-1,1,N);
y = linspace(-1,1,N);
z = linspace(-1,1,N);
[X Y Z] = meshgrid(x,y,z);
dx = x(2)-x(1);

% to return the full matrix, set primed and unprimed coordinates to be the same
Xp = X;
Yp = Y;
Zp = Z;

% build the matrix
G = momGmatrix3D(dx,k,X,Y,Z,Xp,Yp,Zp);

% plot
figure(1),clf
imagesc(real(G))
myplot('MoM Green''s Function Matrix - 3D','Point index (x'',y'',z'')','Point index (x,y,z)')
colorbar
mycbar('real(G)',5e-2*[-1 1])

figure(2),clf
imagesc(imag(G))
myplot('MoM Green''s Function Matrix - 3D','Point index (x'',y'',z'')','Point index (x,y,z)')
colorbar
mycbar('imag(G)',5e-2*[-1 1])

figure(3),clf
imagesc(abs(G))
myplot('MoM Green''s Function Matrix - 3D','Point index (x'',y'',z'')','Point index (x,y,z)')
colorbar
mycbar('abs(G)',5*[0 1e-2])



%% momGmatrixDyadic

k = 2*pi;

% create a dense square grid of points
N = 10;
x = linspace(-1,1,N);
y = linspace(-1,1,N);
z = linspace(-1,1,N);
[X Y Z] = meshgrid(x,y,z);
dx = x(2)-x(1);

% to return the full matrix, set primed and unprimed coordinates to be the same
Xp = X;
Yp = Y;
Zp = Z;

% build the matrix blocks
[Gxx Gyy Gzz Gxy Gxz Gyz] = momGmatrixDyadic(dx,k,X,Y,Z,Xp,Yp,Zp);

% plotting
type = 4;
switch type
    case 1
        tmp = Gxx; lab = 'Gxx'; cx = 3e-2*[0 1];
    case 2
        tmp = Gyy; lab = 'Gyy'; cx = 3e-2*[0 1];
    case 3
        tmp = Gzz; lab = 'Gzz'; cx = 3e-2*[0 1];
    case 4
        tmp = Gxy; lab = 'Gxy'; cx = 1e-2*[0 1];
    case 5
        tmp = Gxz; lab = 'Gxz'; cx = 1e-2*[0 1];
    case 6
        tmp = Gyz; lab = 'Gyz'; cx = 1e-2*[0 1];
end

figure(3),clf
imagesc(abs(tmp))
myplot(['MoM Dyadic Green''s Function Matrix Block ' lab],'Point index (x'',y'',z'')','Point index (x,y,z)')
colorbar
mycbar('abs(G)',cx)









