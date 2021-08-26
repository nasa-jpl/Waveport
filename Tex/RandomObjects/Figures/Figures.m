
direc = '/Users/mshaynes/Desktop/Book/Waveport/Tex/RandomObjects/Figures/';

%% rough1

rng(1)
Nx = 4001;
mid = (Nx+1)/2;
Lx = 30;
rmsh = 2;
lc = 0.2;
x = linspace(-Lx/2,Lx/2,Nx);

surfex = zeros(Nx,2);
cor = zeros(Nx,2);
cor_ana = zeros(Nx,2);
dev = zeros(Nx,2);
dev_ana = zeros(Nx,2);

for m = 1:2,
    
    if m == 1
        type = 'norm';
    else
        type = 'exp';
    end
    
    rng(1)
    % generate surface
    z = rough1(Nx,Lx,rmsh,lc,type);
    surfex(:,m) = z;
    
   
    % compute correlatin and variation 
    Nt = 100;
    data_corr = zeros(Nx,Nt);
    data_dev = zeros(Nx,Nt);

    for t=1:Nt,
        z = rough1(Nx,Lx,rmsh,lc,type);
        data_corr(:,t) = conv(z,flipud(z),'same');
        data_dev(:,t) = z - z(mid);
    end


    cor(:,m) = mean((data_corr),2);
    cor(:,m) = cor(:,m)/max(cor(:,m));
    if m == 1
        cor_ana(:,m) = exp(-x.^2/lc^2);
    else
        cor_ana(:,m) = exp(-abs(x)/lc);
    end
       
    dev(:,m) = rms(data_dev,2);
    if m == 1
        dev_ana(:,m) = sqrt(2)*rmsh*sqrt(1-exp(-x.^2/lc^2));
    else
        dev_ana(:,m) = sqrt(2)*rmsh*sqrt(1-exp(-abs(x)/lc));
    end
end

lim = 10/lc;

h1 = figure(1),clf
plot(x/lc,surfex(:,1)/rmsh)
myplot({'1D Gaussian Random Signal';'Gaussian Correlation Function'},'x/l_c','f(x)/\sigma')
ylim(3*[-1 1])
xlim(lim*[-1 1])


h2 = figure(2),clf
plot(x/lc,surfex(:,2)/rmsh)
myplot({'1D Gaussian Random Signal';'Exponential Correlation Function'},'x/l_c','f(x)/\sigma')
ylim(3*[-1 1])
xlim(lim*[-1 1])


h3 = figure(3),clf,hold all
plot(x/lc,cor_ana(:,1))
plot(x/lc,cor(:,1))
plot(x/lc,cor_ana(:,2))
plot(x/lc,cor(:,2))
hold off
myplot('1D Gaussian Signal Correlation','\Deltax/l_c','C(\Deltax)')
ylim([-0.1 1.1])
legend('Analytic, C_1','Numeric, C_1','Analytic, C_2','Numeric, C_2',...
    'location','northeast')
xlim(lim*[-1 1])


h4 = figure(4),clf,hold all
plot(x/lc,dev_ana(:,1)/rmsh)
plot(x/lc,dev(:,1)/rmsh)
plot(x/lc,dev_ana(:,2)/rmsh)
plot(x/lc,dev(:,2)/rmsh)
hold off
myplot('1D Gaussian Signal RMS Deviation','\Deltax/l_c','RMS Deviation v/\sigma')
ylim([0 1.8])
legend('Analytic, C_1','Numeric, C_1','Analytic, C_2','Numeric, C_2',...
    'location','southeast')
xlim(lim*[-1 1])


if 0
    saveimage(h1,[direc 'gaussiansurf1'],'epsc');
	saveimage(h2,[direc 'gaussiansurf2'],'epsc');
    saveimage(h3,[direc 'gaussiancorr'],'epsc');
	saveimage(h4,[direc 'gaussiandev'],'epsc');

end


%% rough2


Nx = 201;
Ny = 201;
Lx = 40;
Ly = 40;
x = linspace(-Lx/2,Lx/2,Nx);
y = linspace(-Ly/2,Ly/2,Ny);
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
rmsh = 2;
lcx = 1;
lcy = 1;

rng(1);
z1 = rough2(Nx,Ny,Lx,Ly,rmsh,lcx,lcy,'norm');
rng(1);
z2 = rough2(Nx,Ny,Lx,Ly,rmsh,lcx,lcy,'exp');

h1=figure(1),clf
imagesc(x/lcx,y/lcy,z1'/rmsh),colorbar
mycbar('f(x,y)/\sigma',[-3 3])
myplot({'2D Gaussian Rough Surface';'Gaussian Correlation Function'},'x/lcx','y/lcy')
grid off

h2=figure(2),clf
imagesc(x/lcx,y/lcy,z2'/rmsh),colorbar
mycbar('f(x,y)/\sigma',[-3 3])
myplot({'2D Gaussian Rough Surface';'Exponential Correlation Function'},'x/lcx','y/lcy')
grid off


if 0
    saveimage(h1,[direc 'rough2norm'],'epsc');
	saveimage(h2,[direc 'rough2exp'],'epsc');
end



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

h1 = figure(1),clf
slice(X/lcx,Y/lcy,Z/lcz,z1/rmsh,0,0,0),colorbar
mycbar('f(x,y,z)/\sigma',[-3 3])
myplot({'3D Gaussian Random Volume';'Gaussian Correlation Function'},'x/lcx','y/lcy','z/lcz')
shading flat

h2 = figure(2),clf
slice(X/lcx,Y/lcy,Z/lcz,z2/rmsh,0,0,0),colorbar
mycbar('f(x,y,z)/\sigma',[-3 3])
myplot({'3D Gaussian Random Volume';'Exponential Correlation Function'},'x/lcx','y/lcy','z/lcz')
shading flat



if 0
    saveimage(h1,[direc 'gauss3Dnorm'],'epsc');
	saveimage(h2,[direc 'gauss3Dexp'],'epsc');
end



%% bicont

% input parameters
fv = 0.3;       % volume fraction
lave = 5e-3;    % average length scale of heterogeneity
b = 10;         % shape parameter

% evaluation points
Nx = 251;
x = linspace(-0.05,0.05,Nx);
y = x;
z = 0;
[X Y Z] = meshgrid(x,y,z);

% set seed and compute indicator function and fluctuating field
rng(1);
[Theta] = bicont(X,Y,Z,lave,b,fv);

% evaluation points
Nx = 251;
xx = linspace(-0.05,0.05,Nx);
h1 = figure(1),clf,hold all
for cut=1:3,
    switch cut
        case 1
            x = xx;    
            y = xx;
            %Z = z;
            [X Y] = meshgrid(x,y);
            Z = zeros(size(X));
        case 2
            x = xx;    
           % Y = 0;
            z = xx;
            [X Z] = meshgrid(x,z);
            Y = zeros(size(x));
        case 3
            %X = 0;    
            y = xx;
            z = xx;
            [Y Z] = meshgrid(y,z);
            X = zeros(size(Y));
    end
    rng(1);
    [Theta] = bicont(X,Y,Z,lave,b,fv);
    surf(X*1000,Y*1000,Z*1000,Theta)
end

colormap(gray)
shading flat
axis equal
view([1 0.6 0.4])
myplot('Bicontinuous Random Media','x (mm)','y (mm)','z (mm)')
xticks([-50:25:50])
yticks([-50:25:50])
zticks([-50:25:50])



if 0
    saveimage(h1,[direc 'bicont'],'epsc');
end

%% Fractal 1 - RMS Variation

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
    [z] = fractal1(Nx,Lx,H(m),vo(m),dxo,'var');
    surf_examples(:,m) = z;

    surf_rms_analytic(:,m) = vo(m)*(abs(x)/dxo).^H(m);

    Nt = 500;
    data = zeros(Nx,Nt);
    for t=1:Nt,
        [z] = fractal1(Nx,Lx,H(m),vo(m),dxo,'var');
        data(:,t) = z(:);
    end
    surf_rms_statistics(:,m) = rms(data - data(1,:),2);
   
end

h1 = figure(1),clf,hold all
plot(x,surf_examples)
hold off
myplot('1D Fractal Surface Profile','x','Height')
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
myplot({'Fractal Surface RMS Deviation Relative to x = 0';['Number of Trials = ' num2str(Nt)]},'\Deltax','v(\Deltax)')
str1 = ['Surface 1: Analytic'];
str1b = ['Surface 1: Numeric'];
str2 = ['Surface 2: Analytic'];
str2b = ['Surface 2: Numeric'];
str3 = ['Surface 3: Analytic'];
str3b = ['Surface 3: Numeric'];
leg = legend(str1,str1b,str2,str2b,str3,str3b,'location','northwest')




if 0
    saveimage(h1,[direc 'fracsurf'],'epsc');
	saveimage(h2,[direc 'fracdiffvrms'],'epsc');
end



%% Ocean omnidirectional spectrum

u10 = 3;                % wind speed m/s
g = 9.81;               % gravitational acceleration, (m/s^2)
alpha = 0.0081;         % Phillips constant
fp = 0.13*g/u10;      % PM peak freqeucy (Hz)
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


if 0
    saveimage(h1,[direc 'pmspectrum'],'epsc');
end




%%  Ocean PSD

Lx = 100; % domain size in x (m)
Ly = 100; % domain size in y (m)
x = linspace(-Lx/2,Lx/2,400);
y = linspace(-Ly/2,Ly/2,400);
dx = x(2)-x(1);  
dy = y(2)-y(1);
Nx = length(x);
Ny = length(y);
[X Y] = meshgrid(x,y);
wind_direction_deg = 35;     % wind direction (meteorological convention)
[W Kx Ky] = psdOcean(Nx,Ny,Lx,Ly,u10,wind_direction_deg);

h2 = figure(2); clf,
W2 = W/max(max(abs(W)));
imagesc(fftshift(Kx(1,:)),fftshift(Ky(:,1)),fftshift(W2))
myplot2({'Ocean PSD (Normalized)';['Wind Direction ' ...
    num2str(wind_direction_deg) '^o']},'k_x (1/m)','k_y (1/m)')
set(gca,'ydir','normal')
axis square
axis([-2 2 -2 2])
mycbar([],[0 1])


if 0
    saveimage(h2,[direc 'oceanPSD'],'epsc');
end




%% evolveOcean

time = linspace(0,5,50);  % time in seconds
series = zeros([size(W), length(time)]);  % storage for each frame
rng(1)
gam = (1/sqrt(2))*(randn(size(W)) + 1i*randn(size(W)));  % single instance of ocean spectrum

for n=1:length(time),
    
   At = evolveOcean(W,Kx,Ky,gam,time(n));
   At = At/std(At(:))*sigma; % rescale the rms to match the analytical
   
   % plot as running image
   figure(1),clf,imagesc(x,y,At),colorbar
   myplot2({'Wave Heights (m)';['Time: ' num2str(time(n),2) ' (s)']},'x (m)','y (m)')
   axis square
   set(gca,'ydir','normal')
   mycbar('m',[-0.2 0.2])
    
   % save
   series(:,:,n) = At;
   pause(0.01)
   
end



%% plot one ocean surface

h4 = figure(4);clf,
imagesc(x,y,series(:,:,1));
myplot({'Ocean Surface Height';['U10 = ' num2str(u10) ' m/s, Wind Direction = ' num2str(wind_direction_deg) ' deg'];' '},'x (m)','y (m)');
axis square
set(gca,'ydir','normal')
mycbar('m',[-0.2 0.2])
xticks([-50:10:50])
grid off

if 0
    saveimage(h4,[direc 'oceansurfaceheight'],'epsc');
end


%% Gaussian Random Particles

% plot the correlation function

gamma = linspace(0,180,500);
Gamma_degs = [5 15 30 45 60 90 135 180];
[g Gd] = ndgrid(gamma, Gamma_degs);
lc = 2*sind(Gd/2);
C = exp(-2./lc.^2.*sind(g/2).^2);

h1=figure(1),clf,
plot(gamma,C);
hold on
plot([0 180],exp(-1/2)*[1 1],'k--')
hold off
myplot(' ','\gamma (deg)','Correlation')
title('Angular Correlation Function C_s(\gamma)','interpreter','tex')
xlim([0 180])
leg = legend(num2str(Gamma_degs'),'location','eastoutside');
myleg(leg,'\Gamma (deg)')

if 0
    direc = '/Users/mshaynes/Desktop/Book/Chap5/Figures/GaussianParticles/';
    saveimage(h1,[direc 'csgamma'],'epsc');
end



%% plot a few particles

rmshs = [0.03 0.1 0.2];
Gamma_degs = [5 10 20];
GP = struct([]);

for nn=1:length(rmshs),
for gg=1:length(Gamma_degs),
    mydisp(gg,length(Gamma_degs))

    % sphere parameters
    a = 1;          % radius (unitless)
    rmsh = rmshs(nn);
    Gamma_deg = Gamma_degs(gg);
    rng(1);          % set the random seed, or delete

    % theta/phi spacing
    Nt = 200;       % number of theta samples
    Np = 2*Nt;      % number of phi samples = 2x theta samples
    theta = linspace(0,pi,Nt);
    phi = linspace(0,2*pi,Np);
    [Th Phi] = meshgrid(theta,phi);
    [R] = gaussianRandomParticle(Th,Phi,a,rmsh,Gamma_deg);
    [X Y Z] = sph2cart(R,Th,Phi);
    GP(nn,gg).X = X;
    GP(nn,gg).Y = Y;
    GP(nn,gg).Z = Z;
    GP(nn,gg).R = R;

end
end


%% generate figures

cnt = 1;
for gg=1:length(Gamma_degs),
for nn=1:length(rmshs),

    % plot the surface
    h2 = figure(1);,clf

    surf(GP(nn,gg).X,GP(nn,gg).Y,GP(nn,gg).Z,GP(nn,gg).R);
    shading flat
    axis equal
    axis(1.7*a*[-1 1 -1 1 -1 1])
    myplot({'Gaussian Random Particle';...
        ['a = ' num2str(a) ...
        ',   \sigma = ' num2str(rmshs(nn)) ...
        ',   \Gamma = ' num2str(Gamma_degs(gg)) ' deg']}...
        ,'x','y','z')

    pause(0.01)
    if 0
        direc = '/Users/mshaynes/Desktop/Book/Chap5/Figures/GaussianParticles/';
        saveimage(h2,[direc 'part' num2str(cnt)],'jpeg');
    end
    cnt = cnt + 1;


end
end





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
loglog(abs(kx),V)
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
loglog(lamx/1000,V)
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


if 0
    saveimage(h1,[direc 'ionoPSD1'],'epsc');
	saveimage(h2,[direc 'ionoPSD2'],'epsc');
    saveimage(h3,[direc 'ionoProfile'],'epsc');
end





