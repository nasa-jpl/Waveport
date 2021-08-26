 



%% Gaussian quadrature on [-1 1], legpts

figure(1),clf,hold all
all = [];
N = 2:15;
for n=N,
    [xj wj] = legpts(n);
    scatter(xj,wj,40,colors(n),'filled')    
end
hold off
myplot('Nodes and Weights of Gaussian Quadrature','Notes, x_j','Weights, w_j')
leg = legend(num2str(N'),'location','eastoutside')
myleg(leg,'N')
ylim([0 1.1])


%% Spherical nodes of quadrature

L = 6;
I = 2*L+1;
J = L+1;

[xj wj] = legpts(J);
theta_j = acos(xj);
phi_i = (0:(I-1))*2*pi/I;

[P T] = meshgrid(phi_i,theta_j);

figure(2),clf
hold on
plot(180*[-1 1],90*[1 1],'k--')
plot([0 0 ],180*[0 1],'k--')
scatter(180/pi*wrap(P(:)),180/pi*T(:),50,'filled')
hold off
myplot({'Gaussian Quadratue Sampling on Sphere';['L = ' num2str(L)]},'\phi_i (deg)','\theta_j (deg)')
set(gca,'ydir','reverse')
ylim(180/pi*[0 pi])
xlim(180/pi*pi*[-1 1])



%% Translation operator, TLth

Nth = 50;
Nphi = 2*Nth;
th = linspace(0,pi,Nth);
phi = linspace(-pi,pi,Nphi);
[Phi Th] = meshgrid(phi,th);
Ntot = length(Phi(:));

[kxhat kyhat kzhat] = sph2cart(ones(size(Phi)),Th,Phi);
Xhat = [1 0 0];

kX = 50;
costheta = Xhat(1)*kxhat + Xhat(2)*kyhat + Xhat(3)*kzhat;

L = 12;
tl = TLth(L,kX,costheta);

figure(1),clf
imagesc(180/pi*phi,180/pi*th,abs(tl))
myplot({'FMM Translation Operator, TL ';['L = ' num2str(L), ', kX = ' num2str(kX)]},'\phi (deg)','\theta (deg)')
set(gca,'ydir','reverse')
ylim(180/pi*[0 pi])
xlim(180/pi*pi*[-1 1])




%% TL interpolator, interpTLprep

L = 4;
s = 5;
p = 3;
M = s*L;

dth = 2*pi/(2*M+1);
thetao = p*dth;
N = M - L;

r = 5;
k = 2*pi;

[tlsamp, dth, thetao, M, theta_samp] = interpTLprep(L,k,r,s,p);

theta = linspace(0,pi,8*M);
tlinterp = interpTL(theta,L,p,tlsamp,dth,thetao,M);


xlimits = [-0.5 3.75];
ylimits = max(abs(tlsamp))*1.2*[-1 1];
tlsamp = TLth(L,k*r,cos(theta_samp));
h1 = figure(1);,clf,hold all

plot(theta,real(tlinterp),'linewidth',2)
plot(theta,imag(tlinterp),'linewidth',2)
scatter(theta_samp,real(tlsamp),'o','filled')
scatter(theta_samp,imag(tlsamp),'o','filled')

h = plot(xlimits,[0 0],'k','linewidth',2);

plot([0 0],ylimits,'k--','linewidth',2)
plot([pi pi],ylimits,'k--','linewidth',2)
myplot([],'\theta (radians)','Re/Im T_L(\theta)')
legend('Re T_L(\theta) - Interp','Im T_L(\theta) - Interp',...
    'Re T_L(m\Delta\theta) - Sampled','Im T_L(m\Delta\theta) - Sampled')
xlim(xlimits)
ylim(ylimits)
uistack(h,'down',4)


%% sst and isst

L = 30;
I = 2*L+1;
J = L+1;

[muj wj] = legpts(J);
theta_j = acos(muj);
phi_i = (0:(I-1))*2*pi/I;
[T P] = meshgrid(theta_j,phi_i);

% spherical harmonics
ylm = sphericalY(L,T,P,'mono');
tot = L^2 + 2*L + 1;

% random complex expansion coefficients
rng(1)
alm = rand(tot,1) + 1i*rand(tot,1);

% spherical function
f = reshape(ylm*alm,size(T));

% sst
flm = sst(f,L,muj,wj);
mean([alm-flm])

% isst
f2 = isst(flm,L,muj);
mean([f(:)-f2(:)])



%% sst and isst with more points than the number of harmonics (interpolation)

L = 6;
L2 = 10;
I = 2*L2+1;
J = L2+1;
[muj wj] = legpts(J);
theta_j = acos(muj);
phi_i = (0:(I-1))*2*pi/I;
[T P] = meshgrid(theta_j,phi_i);

% spherical harmonics
ylm = sphericalY(L,T,P,'mono');
tot = L^2 + 2*L + 1;

% random complex expansion coefficients
rng(1)
alm = rand(tot,1) + 1i*rand(tot,1);

% spherical function
f = reshape(ylm*alm,size(T));

% sst 
flm = sst(f,L,muj,wj);
% isst
f3 = isst(flm,L,muj);

mean(f3(:)-f(:))

figure(1),clf
imagesc(real(f2));
figure(2),clf
imagesc(real(f3));



%% ssfilt

% field with lower number of harmonics
L = 4;
I = 2*L + 1;
J = L + 1;
phi = 2*pi*(0:(I-1))/I;
[muj, wj] = legpts(J);
theta = acos(muj);
[T1 P1] = meshgrid(theta,phi);
tot1 = L^2 + 2*L + 1;

ylm = sphericalY(L,T1,P1,'mono');
rng(8)
flm = rand(tot1,1) + 1i*rand(tot1,1);
f = reshape(ylm*flm,size(T1));

% field with higher nubmer of harmonics
K = 20;
P = 2*K + 1;
Q = K + 1;
phi2 = 2*pi*(0:(P-1))/P;
[muk, wk] = legpts(Q);
theta2 = acos(muk);
[T2 P2] = meshgrid(theta2,phi2);
tot2 = K^2 + 2*K + 1;

% perform interpolation or filtering
% 0 to interpolate L field to K field
% 1 to filter K field to L field
filter_or_interp = 1;
if filter_or_interp
    ylm = sphericalY(K,T2,P2,'mono');
    rng(3)
    flm = rand(tot2,1) + 1i*rand(tot2,1);
    Lmax = 4;
    cutoff = Lmax^2 + 2*Lmax + 1;
    flm((cutoff+1):end) = 0;
    f = reshape(ylm*flm,size(T2));
end

% function determines whether to filter or interpolate based on size(f),L,K
f2 = ssfilt(f,L,muj,wj,K,muk,wk);

% plot
if ~filter_or_interp
    
    h1 = figure(1);clf,imagesc(180/pi*phi,180/pi*theta,real(f)'),colorbar,caxis([-1 3])
    myplot(['L = ' num2str(L)],'\phi (deg)','\theta (deg)')
    grid off
    
    h2 = figure(2);clf,imagesc(180/pi*phi2,180/pi*theta2,real(f2)'),colorbar,caxis([-1 3])
    myplot(['K = ' num2str(K)],'\phi (deg)','\theta (deg)')
    grid off
    
else
    
    h34 = figure(3);clf,imagesc(180/pi*phi,180/pi*theta,real(f)'),colorbar,caxis([-2 2])
    myplot(['K = ' num2str(K)],'\phi (deg)','\theta (deg)')
    grid off
    
    h4 = figure(4);clf,imagesc(180/pi*phi2,180/pi*theta2,real(f2)'),colorbar,caxis([-2 2])
    myplot(['L = ' num2str(L)],'\phi (deg)','\theta (deg)')
    grid off
end


%% fssfilt

% field with lower number of harmonics
L = 10;
I = 2*L + 1;
J = L + 1;
phi = 2*pi*(0:(I-1))/I;
[muj, wj] = legpts(J);
theta = acos(muj);
[T1 P1] = meshgrid(theta,phi);
tot1 = L^2 + 2*L + 1;

ylm = sphericalY(L,T1,P1,'mono');
rng(8)
flm = rand(tot1,1) + 1i*rand(tot1,1);
f = reshape(ylm*flm,size(T1));

% field with higher nubmer of harmonics
K = 30;
P = 2*K + 1;
Q = K + 1;
phi2 = 2*pi*(0:(P-1))/P;
[muk, wk] = legpts(Q);
theta2 = acos(muk);
[T2 P2] = meshgrid(theta2,phi2);
tot2 = K^2 + 2*K + 1;

% perform interpolation or filtering
% 0 to interpolate L field to K field
% 1 to filter K field to L field
filter_or_interp = 1;
if filter_or_interp
    ylm = sphericalY(K,T2,P2,'mono');
    rng(3)
    flm = rand(tot2,1) + 1i*rand(tot2,1);
    Lmax = 4;
    cutoff = Lmax^2 + 2*Lmax + 1;
    flm((cutoff+1):end) = 0;
    f = reshape(ylm*flm,size(T2));
end

% function determines whether to filter or interpolate based on size(f),L,K
tic
f2 = ssfilt(f,L,muj,wj,K,muk,wk);
toc
tic
f3 = fssfilt(f,L,muj,wj,K,muk,wk);
toc
mean(f2(:)-f3(:))

% Try with lower harmonics than sampling
f2 = ssfilt(f,L-1,muj,wj,K-1,muk,wk);
f3 = fssfilt(f,L-1,muj,wj,K-1,muk,wk);
mean(f2(:)-f3(:))


%% vst 

L = 10;
I = 2*L + 1;
J = L + 1;
phi = 2*pi*(0:(I-1))/I;
[muj, wj] = legpts(J);
theta = acos(muj);
[T P] = meshgrid(theta,phi);
tot = L^2 + 2*L;

% partially normalized B and C
[Bth, Bphi, Cth, Cphi] = BC(L,T,P);

rng(1);
blm1 = rand(tot,1) + 1i*rand(tot,1);
clm1 = rand(tot,1) + 1i*rand(tot,1);
[Fth Fphi] = BCmult(Bth,Bphi,Cth,Cphi,blm1,clm1);
Fth = reshape(Fth,size(T));
Fphi = reshape(Fphi,size(T));

% vst to compute blm and clm
[blm, clm] = vst(Fth,Fphi,L,muj,wj);
mean(abs(blm1-blm))
mean(abs(clm1-clm))

% vst with precomputation 
plm = Plm(L,muj); 
dplm = Plmp(L,muj);
[blm, clm] = vst(Fth,Fphi,L,muj,wj,[],plm,dplm);
mean(abs(blm1-blm))
mean(abs(clm1-clm))

% fully normalized B and C
[Bth, Bphi, Cth, Cphi] = BC(L,T,P,'norm');
[Fth Fphi] = BCmult(Bth,Bphi,Cth,Cphi,blm1,clm1);
Fth = reshape(Fth,size(T));
Fphi = reshape(Fphi,size(T));

[blm, clm] = vst(Fth,Fphi,L,muj,wj,'norm');
mean(abs(blm1-blm))
mean(abs(clm1-clm))

% fully normalized vst with precomputation 
plm = Plm(L,muj); 
dplm = Plmp(L,muj);
[blm, clm] = vst(Fth,Fphi,L,muj,wj,'norm',plm,dplm);
mean(abs(blm1-blm))
mean(abs(clm1-clm))



%% ivst, L matches I and J sampling

L = 10;
I = 2*L + 1;
J = L + 1;
phi = 2*pi*(0:(I-1))/I;
[muj, wj] = legpts(J);
theta = acos(muj);
[T P] = meshgrid(theta,phi);
tot = L^2 + 2*L;

% partially normalized B and C
[Bth, Bphi, Cth, Cphi] = BC(L,T,P);

rng(1);
blm = rand(tot,1) + 1i*rand(tot,1);
clm = rand(tot,1) + 1i*rand(tot,1);
[Fth Fphi] = BCmult(Bth,Bphi,Cth,Cphi,blm,clm);
Fth = reshape(Fth,size(T));
Fphi = reshape(Fphi,size(T));

% partially normalized
[Fth2, Fphi2] = ivst(blm,clm,L,muj);
mean(abs(Fth(:)-Fth2(:)))
mean(abs(Fphi(:)-Fphi2(:)))

% test with precomputation, partially normalized
plm = Plm(L,muj); 
dplm = Plmp(L,muj);
[Fth2, Fphi2] = ivst(blm,clm,L,muj,[],plm,dplm);
mean(abs(Fth(:)-Fth2(:)))
mean(abs(Fphi(:)-Fphi2(:)))

% fully normalized B and C
[Bth, Bphi, Cth, Cphi] = BC(L,T,P,'norm');
[Fth Fphi] = BCmult(Bth,Bphi,Cth,Cphi,blm,clm);
Fth = reshape(Fth,size(T));
Fphi = reshape(Fphi,size(T));

[Fth2, Fphi2] = ivst(blm,clm,L,muj,'norm');
mean(abs(Fth(:)-Fth2(:)))
mean(abs(Fphi(:)-Fphi2(:)))

% fully normalized B and C
plm = Plm(L,muj); 
dplm = Plmp(L,muj);
[Fth2, Fphi2] = ivst(blm,clm,L,muj,'norm',plm,dplm);
mean(abs(Fth(:)-Fth2(:)))
mean(abs(Fphi(:)-Fphi2(:)))


%% ivst, I and J are larger than L

% sampling at L2 harmonics
L2 = 15;
I = 2*L2 + 1;
J = L2 + 1;
phi = 2*pi*(0:(I-1))/I;
[muj, wj] = legpts(J);
theta = acos(muj);
[T P] = meshgrid(theta,phi);

% fewer than harmonics than the sampling
L = 10;
tot = L^2 + 2*L;

% partially normalized B and C
[Bth, Bphi, Cth, Cphi] = BC(L,T,P);

rng(1);
blm = rand(tot,1) + 1i*rand(tot,1);
clm = rand(tot,1) + 1i*rand(tot,1);

[Fth Fphi] = BCmult(Bth,Bphi,Cth,Cphi,blm,clm);
Fth = reshape(Fth,size(T));
Fphi = reshape(Fphi,size(T));

% I and J are sampled for L2 harmonics, partially normalized B and C
[Fth2,Fphi2] = ivst(blm,clm,L,muj);
mean(abs(Fth(:)-Fth2(:)))
mean(abs(Fphi(:)-Fphi2(:)))

% I and J are sampled for L2 harmonics, fully normalized B and C
[Bth, Bphi, Cth, Cphi] = BC(L,T,P,'norm');
[Fth Fphi] = BCmult(Bth,Bphi,Cth,Cphi,blm,clm);
Fth = reshape(Fth,size(T));
Fphi = reshape(Fphi,size(T));

[Fth2,Fphi2] = ivst(blm,clm,L,muj,'norm');
mean(abs(Fth(:)-Fth2(:)))
mean(abs(Fphi(:)-Fphi2(:)))



%% vsfilt

% field with lower number of harmonics
L = 12;
I = 2*L + 1;
J = L + 1;
phi = 2*pi*(0:(I-1))/I;
[muj, wj] = legpts(J);
theta = acos(muj);
[T1 P1] = meshgrid(theta,phi);
tot1 = L^2 + 2*L;

[Bth, Bphi, Cth, Cphi] = BC(L,T1,P1);
rng(3)
blm = rand(tot1,1)-0.5 + 1i*(rand(tot1,1)-0.5);
clm = rand(tot1,1)-0.5 + 1i*(rand(tot1,1)-0.5);
[Fth Fphi] = BCmult(Bth, Bphi, Cth, Cphi,blm,clm);
Fth = reshape(Fth,size(T1));
Fphi = reshape(Fphi,size(T1));


% field with higher nubmer of harmonics
K = 30;
P = 2*K + 1;
Q = K + 1;
phi2 = 2*pi*(0:(P-1))/P;
[muk, wk] = legpts(Q);
theta2 = acos(muk);
[T2 P2] = meshgrid(theta2,phi2);
tot2 = K^2 + 2*K;

% perform interpolation or filtering
% 0 to interpolate L field to K field
% 1 to filter K field to L field
filter_or_interp = 1;
if filter_or_interp
    [Bth, Bphi, Cth, Cphi] = BC(K,T2,P2);
    rng(3)
    blm = rand(tot2,1) + 1i*rand(tot2,1);
    clm = rand(tot2,1) + 1i*rand(tot2,1);
    Lmax = 4;
    cutoff = Lmax^2 + 2*Lmax;
    blm((cutoff+1):end) = 0;
    clm((cutoff+1):end) = 0;
    
    [Fth Fphi] = BCmult(Bth, Bphi, Cth, Cphi,blm,clm);
    Fth = reshape(Fth,size(T2));
    Fphi = reshape(Fphi,size(T2));
end

% function determines whether to filter or interpolate based on size(f),L,K
[Fth2 Fphi2] = vsfilt(Fth,Fphi,L,muj,wj,K,muk,wk);

% plot
if ~filter_or_interp
    
    h1 = figure(1);clf,imagesc(180/pi*phi,180/pi*theta,imag(Fth)'),colorbar,caxis(5*[-1 1])
    myplot(['L = ' num2str(L)],'\phi (deg)','\theta (deg)')
    grid off
    
    h2 = figure(2);clf,imagesc(180/pi*phi,180/pi*theta,imag(Fphi)'),colorbar,caxis(5*[-1 1])
    myplot(['L = ' num2str(L)],'\phi (deg)','\theta (deg)')
    grid off
    
    h3 = figure(3);clf,imagesc(180/pi*phi2,180/pi*theta2,imag(Fth2)'),colorbar,caxis(5*[-1 1])
    myplot(['K = ' num2str(K)],'\phi (deg)','\theta (deg)')
    grid off
    
    h4 = figure(4);clf,imagesc(180/pi*phi2,180/pi*theta2,imag(Fphi2)'),colorbar,caxis(5*[-1 1])
    myplot(['K = ' num2str(K)],'\phi (deg)','\theta (deg)')
    grid off
        
else
    
    h1 = figure(1);clf,imagesc(180/pi*phi,180/pi*theta,real(Fth)'),colorbar,caxis(5*[-1 1])
    myplot(['K = ' num2str(K)],'\phi (deg)','\theta (deg)')
    grid off
    
    h2 = figure(2);clf,imagesc(180/pi*phi,180/pi*theta,real(Fphi)'),colorbar,caxis(5*[-1 1])
    myplot(['K = ' num2str(K)],'\phi (deg)','\theta (deg)')
    grid off
    
    h3 = figure(3);clf,imagesc(180/pi*phi2,180/pi*theta2,real(Fth2)'),colorbar,caxis(5*[-1 1])
    myplot(['L = ' num2str(L)],'\phi (deg)','\theta (deg)')
    grid off
    
    h4 = figure(4);clf,imagesc(180/pi*phi2,180/pi*theta2,real(Fphi2)'),colorbar,caxis(5*[-1 1])
    myplot(['L = ' num2str(L)],'\phi (deg)','\theta (deg)')
    grid off
    
end


%% vsfilt: fewer harmonics than the sampling, using settings above

% function determines whether to filter or interpolate based on size(f),L,K
[Fth2 Fphi2] = vsfilt(Fth,Fphi,L-1,muj,wj,K-1,muk,wk);


%% fvsfilt1, fvsfilt2: fast vector spherical filter - methods 1 and 2, using settings above

% using the setup above
if filter_or_interp
    str = 'filter';
else
    str = 'interp';
end

% Method 1
[Am, Bm] = fvsfilt1AmBm(L,K,str);
[Fth3 Fphi3] = fvsfilt1(Fth,Fphi,L,K,Am,Bm);
% compare to vsfilt
mean(abs(Fth2(:) - Fth3(:)))
mean(abs(Fphi2(:) - Fphi3(:)))

% Method 2
[Fth4 Fphi4] = fvsfilt2(Fth,Fphi,L,muj,wj,K,muk,wk);
% compare to vsfilt
mean(abs(Fth2(:) - Fth4(:)))
mean(abs(Fphi2(:) - Fphi4(:)))

%%% check the timing between all three
tic
[Fth2 Fphi2] = vsfilt(Fth,Fphi,L,muj,wj,K,muk,wk);
toc
tic
[Fth3 Fphi3] = fvsfilt1(Fth,Fphi,L,K,Am,Bm);
toc
tic
[Fth4 Fphi4] = fvsfilt2(Fth,Fphi,L,muj,wj,K,muk,wk);
toc




%% S-matrix to T-matrix transformations with quadrature

L = 10;
I = 2*L + 1;
J = L + 1;
k = 1;
N = L^2 + 2*L;


% ficticious T-matrix
rng(1);
sz = N;
Tmm = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);
Tmn = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);
Tnm = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);
Tnn = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);

% T to S
[Stt, Stp, Spt, Spp] = convert_T_to_S(Tmm,Tmn,Tnm,Tnn,k,I,J,L);

figure(1),clf
imagesc(real(Stt(:,:,1,1)))
imagesc(real(Stp(:,:,1,1)))
imagesc(real(Spt(:,:,1,1)))
imagesc(real(Spp(:,:,1,1)))

%%% S to T 
[Tmm2, Tmn2, Tnm2, Tnn2] = convert_S_to_T(Stt,Stp,Spt,Spp,k,I,J,L);

mean(abs(Tmm(:) - Tmm2(:)))
mean(abs(Tmn(:) - Tmn2(:)))
mean(abs(Tnm(:) - Tnm2(:)))
mean(abs(Tnn(:) - Tnn2(:)))

%% S-matrix integration with quadrature

L = 10;
I = 2*L + 1;
J = L + 1;
k = 1;
N = L^2 + 2*L;

% ficticious T-matrix
rng(1);
sz = N;
Tmm = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);
Tmn = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);
Tnm = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);
Tnn = (rand(sz)-0.5) + 1i*(rand(sz)-0.5);

% T to S
[Stt, Stp, Spt, Spp] = convert_T_to_S(Tmm,Tmn,Tnm,Tnn,k,I,J,L);

% incident field far field coefficients
blm = (rand(N,1)-0.5) + 1i*(rand(N,1)-0.5);
clm = (rand(N,1)-0.5) + 1i*(rand(N,1)-0.5);

% compute incidend field pattern from coefficients
[muj, wj] = legpts(J);
[Gth Gphi] = ivst(blm,clm,L,muj,'norm');

% apply scattering matrix
[Fth, Fphi] = compute_Smatrix_quad(Gth,Gphi,Stt,Stp,Spt,Spp,I,J);

% far-field T-matrix 
[Tbb, Tbc, Tcb, Tcc] = convert_T_to_Tfar(Tmm,Tmn,Tnm,Tnn,L,k);

% far-field scattered field from the far-field T-matrix
blmp = Tbb*blm + Tbc*clm;
clmp = Tcb*blm + Tcc*clm;

% scattered field pattern from coefficents
[Fth2 Fphi2] = ivst(blmp,clmp,L,muj,'norm');

% compare fields from S-matrix intergal and far-field T-matrix
mean(abs(Fth(:) - Fth2(:)))
mean(abs(Fphi(:) - Fphi2(:)))



%% 1D FMM

% number of souces
Nk = 40; 
% source points
rng(5);
xk = sort(2*(rand(Nk,1)-0.5));
% source amplitudes
ak = ones(Nk,1); 

% number of observation points
Nj = 10000;  
% observation points
a = -1;
b = 1;
xj = linspace(a,b,Nj)'; 

% precision 
prec = 1e-5;

% prepatory function
[S] = fmm1prep(xk,xj,prec);

% 1D fast multipole method
tic
[f] = fmm1(ak,S);
toc

% direct matrix solution
[Xj Xk] = ndgrid(xj,xk);
M = 1./(Xj-Xk);
ind = find(M == Inf);
M(ind) =0;

tic
fdir = M*ak;
toc


figure(1),clf,hold all
plot(xj,fdir)
plot(xj,f,':')
scatter(xk,zeros(Nk,1),100,'kx')
hold off
legend('f(x_j) Direct','f(x_j) 1D FMM','Sources')
myplot('Direct and 1D FMM solution','x','f(x_j)')
ylim(1e3*[-1 1])

figure(2),clf,
plot(xj,[fdir-f])
myplot('Residual of Direct and 1D FMM','x','f_d - f_{fmm1}')
ylim(2*prec*[-1 1])







