




%% rough1


Nx = 2001;
Lx = 30;
rmsh = 2;
lc = 1;
x = linspace(-Lx/2,Lx/2,Nx);
dx = Lx/(Nx-1);


z1 = rough1(Nx,Lx,rmsh,lc,'norm');
z2 = rough1(Nx,Lx,rmsh,lc,'exp');

rms(z1)
rms(z2)

figure(1),clf
plot(x/lc,z1/rmsh);
myplot({'1D Rough Surface';'Gaussian Correlatin Function'},'x/lc','z/\sigma')
xlim([x(1) x(end)]/lc)
ylim([-3 3])

figure(2),clf
plot(x/lc,z2/rmsh);
myplot({'1D Rough Surface';'Exponential Correlatin Function'},'x/lc','z/\sigma')
xlim([x(1) x(end)]/lc)
ylim([-3 3])



%% rough2


Nx = 301;
Ny = 301;
Lx = 100;
Ly = 100;
rmsh = 2;
lcx = 2;
lcy = 2;
x = linspace(-Lx/2,Lx/2,Nx);
y = linspace(-Ly/2,Ly/2,Ny);
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);


rng(1);
z1 = rough2(Nx,Ny,Lx,Ly,rmsh,lcx,lcy,'norm');
rng(1);
z2 = rough2(Nx,Ny,Lx,Ly,rmsh,lcx,lcy,'exp');

rms(z1(:))
rms(z2(:))

figure(1),clf
imagesc(x/lcx,y/lcy,z1/rmsh),colorbar
mycbar('f(x,y)/\sigma',[-3 3])
myplot({'2D Gaussian Rough Surface';'Gaussian Correlation Function'},'x/lcx','y/lcy')
set(gca,'ydir','normal')
grid off

figure(2),clf
imagesc(x/lcx,y/lcy,z2/rmsh),colorbar
mycbar('f(x,y)/\sigma',[-3 3])
myplot({'2D Gaussian Rough Surface';'Exponential Correlation Function'},'x/lcx','y/lcy')
set(gca,'ydir','normal')
grid off



%% rough3


Nx = 71;
Ny = 71;
Nz = 71;
Lx = 40;
Ly = 40;
Lz = 40;
rmsh = 2;
lcx = 2;
lcy = 2;
lcz = 2;
x = linspace(-Lx/2,Lx/2,Nx);
y = linspace(-Ly/2,Ly/2,Ny);
z = linspace(-Lz/2,Lz/2,Nz);
[X Y Z] = meshgrid(x,y,z);
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
dz = Lz/(Nz-1);


rng(1);
z1 = rough3(Nx,Ny,Nz,Lx,Ly,Lz,rmsh,lcx,lcy,lcz,'norm');
rng(1);
z2 = rough3(Nx,Ny,Nz,Lx,Ly,Lz,rmsh,lcx,lcy,lcz,'exp');

rms(z1(:))
rms(z2(:))

figure(1),clf
slice(X/lcx,Y/lcy,Z/lcz,z1/rmsh,0,0,0),colorbar
mycbar('f(x,y,z)/\sigma',[-3 3])
myplot({'3D Gaussian Random Volume';'Gaussian Correlation Function'},'x/lcx','y/lcy','z/lcz')
shading flat

figure(2),clf
slice(X/lcx,Y/lcy,Z/lcz,z2/rmsh,0,0,0),colorbar
mycbar('f(x,y,z)/\sigma',[-3 3])
myplot({'3D Gaussian Random Volume';'Exponential Correlation Function'},'x/lcx','y/lcy','z/lcz')
shading flat



%% bicont

% input parameters
fv = 0.3;       % volume fraction
lave = 5e-3;    % average length scale of heterogeneity
b = 10;         % shape parameter

% compute gamma distribution parameters
scale = kave/(b+1);
shape = b+1;

% compute cutting level
alp = erfinv(1 - 2*fv);

% evaluation points
Nx = 251;
x = linspace(-0.05,0.05,Nx);
y = x;
z = 0;
[X Y Z] = meshgrid(x,y,z);

% set seed and compute indicator function and fluctuating field
rng(1);
[Theta S] = bicont(X,Y,Z,lave,b,fv);

% check the statistics
mean(S(:))
var(S(:))

% plot
figure(1),clf,
imagesc(x*1000,y*1000,S),colorbar
colormap(gray)
myplot('Bicontinous Media Fluctuating Field','x (mm)','y (mm)')

figure(2),clf,
imagesc(x*1000,y*1000,Theta)
colormap(gray)
myplot('Bicontinous Media Indicator Function','x (mm)','y (mm)')

% check that the mean of the indicator function equals the volume fraction
mean(Theta(:))






%% Fractal 1 - RMS Deviation - Roll off

h = 25e3;
Nx = 2001;
Lx = 100;
x = linspace(0,Lx,Nx)';
 
surf_examples = zeros(Nx,2);
surf_rms_analytic = zeros(Nx,2);
surf_rms_statistics = zeros(Nx,2);

dxo = 1;
vo = [1 1 1];
H = [0.3 0.5 0.7];

for m = 1:length(H),
   
    rng(4)
    [z] = fractal1(Nx,Lx,H(m),vo(m),dxo,'dev');
    surf_examples(:,m) = z;

    surf_rms_analytic(:,m) = vo(m)*(abs(x)/dxo).^H(m);

    Nt = 500;
    data = zeros(Nx,Nt);
    for t=1:Nt,
        [z] = fractal1(Nx,Lx,H(m),vo(m),dxo,'dev');
        data(:,t) = z(:);
    end
    surf_rms_statistics(:,m) = rms(data - data(1,:),2);
   
end

h1 = figure(1),clf,hold all
plot(x,surf_examples)
hold off
myplot('Cross Track Fractal Surface Profile','x','Height')
str1 = ['Surface 1: H = ' num2str(H(1)) ', v_{o} = ' num2str(vo(1)) ', ' '\Deltax_o = ' num2str(dxo)];
str2 = ['Surface 2: H = ' num2str(H(2)) ', v_{o} = ' num2str(vo(2)) ', ' '\Deltax_o = ' num2str(dxo)];
str3 = ['Surface 3: H = ' num2str(H(3)) ', v_{o} = ' num2str(vo(3)) ', ' '\Deltax_o = ' num2str(dxo)];
leg = legend(str1,str2,str3,'location','north')

h2 = figure(2),clf,hold all
plot(x,surf_rms_analytic(:,1))
plot(x,surf_rms_statistics(:,1))
plot(x,surf_rms_analytic(:,2))
plot(x,surf_rms_statistics(:,2))
plot(x,surf_rms_analytic(:,3))
plot(x,surf_rms_statistics(:,3))
hold off
myplot({'Fractal Surface RMS Deviationi Relative to x = 0';['Number of Trials = ' num2str(Nt)]},'x','v(\Deltax)')
str1 = ['Surface 1: Analytic'];
str1b = ['Surface 1: Numeric'];
str2 = ['Surface 2: Analytic'];
str2b = ['Surface 2: Numeric'];
str3 = ['Surface 3: Analytic'];
str3b = ['Surface 3: Numeric'];
leg = legend(str1,str1b,str2,str2b,str3,str3b,'location','northwest')



%% Ocean omnidirectional spectrum

u10 = 3;                % wind speed m/s
g = 9.81;               % gravitational acceleration, (m/s^2)
alpha = 0.0081;         % Phillips constant
fp = 0.13*g/u10;        % PM peak freqeucy (Hz)
const_f = alpha*g^2*(2*pi)^-4;   % multiplying constant used throughout

% basic PM spectrum
f = linspace(0.01,2,1000);   % frequency (Hz)
F = const_f*(f.^-5).*exp(-(5/4)*(fp./f).^4);  % PM spectrum

% RMS height 
sigma = sqrt(alpha)*g/(2*pi)^2/sqrt(5)/fp^2;

% plot spectrum vs f and lambda
h1 = figure(1),clf,
subplot(211),plot(f,F)
myplot('Pierson - Moskowitz Omnidirection Ocean Spectrum', ...
    'Frequency (Hz)','F(f) (m^2/s)')
axis([0 1.5 0 0.013])
subplot(212),plot(g./(2*pi*f.^2),F)
myplot([],'Wavelength (m)','F(f) (m^2/s)')
axis([0.1 30 0 0.012])


%% psdOcean

Lx = 100; % domain size in x (m)
Ly = 100; % domain size in y (m)
x = linspace(-Lx/2,Lx/2,500);
y = linspace(-Ly/2,Ly/2,500);
dx = x(2)-x(1);  
dy = y(2)-y(1);
Nx = length(x);
Ny = length(y);
[X Y] = meshgrid(x,y);

wind_direction_deg = 35;     % wind direction (meteorological convention)
[W Kx Ky] = psdOcean(Nx,Ny,Lx,Ly,u10,wind_direction_deg);

h10 = figure(10); clf,
W2 = W/max(max(abs(W)));
imagesc(fftshift(Kx(1,:)),fftshift(Ky(:,1)),fftshift(W2))
myplot2({'Ocean Spectrum Magnitude (Normalized)';['Wind Direction ' num2str(wind_direction_deg) '^o']},'k_x (1/m)','k_y (1/m)')
set(gca,'ydir','normal')
axis square
axis([-2 2 -2 2])
mycbar([],[0 1])




%% evolveOcean

time = linspace(0,5,100);  % time in seconds
series = zeros([size(W), length(time)]);  % storage for each frame
rng(1)
gam = (1/sqrt(2))*(randn(size(W)) + 1i*randn(size(W)));  % single instance of ocean spectrum

for n=1:length(time),
    
   At = evolveOcean(W,Kx,Ky,gam,time(n));
   At = At/std(At(:))*sigma;
   % plot as running image
   figure(1),clf,imagesc(x,y,At),colorbar
   myplot2({'Wave Heights (m)';['Time: ' num2str(time(n),2) ' (s)']},'x (m)','y (m)')
   axis square
   set(gca,'ydir','normal')
   mycbar('m',[-0.2 0.2])
    
   % save the series
   series(:,:,n) = At;
   pause(0.01)
   
end


%% Gaussian Random Particles

% double check the paper's log-normal and Gaussian transform
N = 1e5;
sigma = 2;
a = 10;
beta = sqrt(log(sigma^2 + 1));
s = beta*randn(N,1);
r = a/sqrt(1 + sigma^2)*exp(s);
[0 mean(s) beta std(s)]
[a mean(r) a*sigma std(r)]


%% make a particle

% sphere parameters
a = 1;          % radius (unitless)
rmsh = 0.03;  	% log-normal RMS of the surface (unitless)
Gamma_deg = 8;  % angular correlation, in degrees.  Gamma_deg > 2.2 deg at the moment
rng(1);          % set the random seed, or delete

% theta/phi spacing
Nt = 100;       % number of theta samples
Np = 2*Nt;      % number of phi samples = 2x theta samples
theta = linspace(0,pi,Nt);
phi = linspace(0,2*pi,Np);
[Th Phi] = meshgrid(theta,phi);

[R] = gaussianRandomParticle(Th,Phi,a,rmsh,Gamma_deg);

[X Y Z] = sph2cart(R,Th,Phi);

% plot the surface
surf(X,Y,Z,R);
shading flat
axis equal
axis(1.4*a*[-1 1 -1 1 -1 1])
myplot({'Gaussian Random Particle';...
    ['Radius: a = ' num2str(a)]; ...
    ['log-Normal RMS: a\sigma = ' num2str(rmsh)]; ...
    ['Angular correlation: \Gamma = ' num2str(Gamma_deg) ' deg']}...
    ,'x','y','z')
mycbar('Radius',a + 4*[-rmsh rmsh])




%% ionosphere1

% number of sampling points
Nx = 5000;

% length of profile in (m)
Lx = 100e3;

% axis
x = linspace(0,Lx,Nx);
dx = x(2)-x(1);

% avearage and percent RMS (here for Ne)
ave = 1e12;
pctrms = 3;

rng(1)
version = 1;
switch version 
    case 1 % use default values
        [Ne V kx Lo Lb] = ionosphere1(Nx,Lx,ave,pctrms);
    case 2 % user defined inputs
        Lo = 20e3;
        Lb = 2000;
        v1 = 3/2;
        v2 = 5/2;
        [Ne V kx] = ionosphere1(Nx,Lx,ave,pctrms,Lo,Lb,v1,v2);
end

% check rms scaling: subtract mean, compute RMS, divide by mean to get percent)
rms(Ne - mean(Ne))/mean(Ne)


% outerscale and breakscale wavenumbers
ko = 2*pi/Lo;
kb = 2*pi/Lb;

% set plotting limits
ylims = [min(V)/10 max(V)*10];
xlims_m = [0.99*dx 2*Lx];


% plot PSD vs wavenumber
h1 = figure(1),clf
loglog(abs(kx),V,'o')
hold all
plot(ko*[1 1],ylims,'k--')
plot(kb*[1 1],ylims,'k--')
hold off
myplot('2-Parameter PSD of Ionosphere Irregularity','Wavenumber, k = 2\pi/L, (1/m)','PSD')
xlim(2*pi./fliplr(xlims_m))
ylim(ylims)

% plot PSD vs lengthscale
lamx = 2*pi./abs(kx);
h2 = figure(2),clf
loglog(lamx/1000,V,'o')
hold all
plot(Lo*[1 1]/1000,ylims,'k--')
plot(Lb*[1 1]/1000,ylims,'k--')
hold off
myplot('2-Parameter PSD of Ionosphere Irregularity','Length scale, L (km)','PSD')
xlim(xlims_m/1000)
ylim(ylims)

% plot the profile
h3 = figure(3),clf
plot(x/1000,Ne)
ylim(ave*[0.8 1.2])
myplot({'Ionosphere Irregularity';...
    ['Average = ' num2str(ave,2) ' el/m^3,   RMS = ' num2str(pctrms) '%']},...
    'Distance (km)','Electron Density (el/m^3)')




