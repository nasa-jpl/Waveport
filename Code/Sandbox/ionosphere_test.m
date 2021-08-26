

Nx = 5000;
Lx = 100e3;

x = linspace(0,Lx,Nx);
dx = x(2)-x(1);

% set avearage and percent RMS
ave = 1e12;
pctrms = 3;

rng(1)
version = 1;
switch version 
    case 1 % defaults
        [Ne V kx Lo Lb] = ionosphere1(Nx,Lx,ave,pctrms);
    case 2 % user defined inputs
        Lo = 20e3;
        Lb = 2000;
        v1 = 3/2;
        v2 = 5/2;
        [Ne V kx] = ionosphere1(Nx,Lx,ave,pctrms,Lo,Lb,v1,v2);
end

% check rms: subtract mean, compute RMS, divide by mean to get percent)
rms(Ne - mean(Ne))/mean(Ne)


figure(1),clf
plot(x/1000,Ne)
ylim(ave*[0.8 1.2])
myplot({'Ionosphere Irregularity';...
    ['Average = ' num2str(ave,2) ' el/m^3,   RMS = ' num2str(pctrms) '%']},...
    'Distance (km)','Electron Density (el/m^3)')



ko = 2*pi/Lo;
kb = 2*pi/Lb;

ylims = [min(V)/10 max(V)*10];
xlims_m = [0.99*dx 2*Lx];

figure(2),clf
loglog(abs(kx),V)
hold all
plot(ko*[1 1],ylims,'k--')
plot(kb*[1 1],ylims,'k--')
hold off
myplot('2-Parameter PSD of Ionosphere Irregularity','Wavenumber, k = 2\pi/L, (1/m)','Power (dB)')
xlim(2*pi./fliplr(xlims_m))
ylim(ylims)


lamx = 2*pi./abs(kx);
figure(3),clf
loglog(lamx/1000,V)
hold all
plot(Lo*[1 1]/1000,ylims,'k--')
plot(Lb*[1 1]/1000,ylims,'k--')
hold off
myplot('2-Parameter PSD of Ionosphere Irregularity','Length scale, L (km)','Spectral Power (dB)')
xlim(xlims_m/1000)
ylim(ylims)


%% compute correlation function (wrapped)


Ntrials = 2000;
Cs = zeros(Nx,Ntrials);

for n=1:Ntrials,

    [Ne V kx] = ionosphere1(Nx,Lx,pctrms,ave);

    Ne2 = Ne-mean(Ne);
    F = fft(Ne2);
    C = fftshift(ifft(F.*conj(F)));
    C = C./max(C);
    Cs(:,n) = C(:);

end
C = mean(Cs,2);

x2 = (x-mean(x))/1000;
f = exp(-abs(x2)/1);
figure(4),clf,hold all
plot(x2-dx/2/1000,C)
plot(x2,f)
hold off


