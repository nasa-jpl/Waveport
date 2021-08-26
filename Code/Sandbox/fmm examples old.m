 



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
Xhat = unit(Xhat);

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
rng(2);
xk = sort(2*(rand(Nk,1)-0.5));
% source amplitudes
ak = ones(Nk,1); 

% number of observation points
Nj = 500;  
% observation points
xj = linspace(a,b-step,Nj)'; 

% precision 
prec = 1e-5;

% prepatory function
tic 
[S] = fmm1prep(xk,xj,prec);
toc 

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


figure(1),clf,
plot(xj,[fdir f])
myplot('Direct and 1D FMM solution','x','f(x_j)')
ylim(1e3*[-1 1])

figure(2),clf,
plot(xj,[fdir-f])
myplot('Residual of Direct and 1D FMM','x','f_d - f_{fmm1}')
ylim(2*prec*[-1 1])





%%

k1 = 2*pi;
d = 10;
ep = 15;

kd = k1*d;

prec = 1e-5;
aph = log10(1/prec);
L = kd + 1.8*aph^(2/3)*kd^(1/3);
L = round(L);
tot = L^2 + 2*L + 1
tot2 = (2*L+1)*(L+1)


[tmm tnn] = sphereTmatrix(L,d,k1,k1/sqrt(ep),1,1);

plot(10*log10(abs([tmm tnn])))
%%


d = 1:200;
prec = logspace(-1,-12,20);
[D P] = meshgrid(d,prec);

kd = 2*pi*D;
aph = log10(1./P);
L = ceil(kd + 1.8*aph.^(2/3).*kd.^(1/3));

%%

h7 = figure(7),clf,contour(d,-log10(prec),L,20),colorbar
set(gca,'ydir','normal')
myplot('Translation Operator Precision vs Truncation Degree','Group Dimension (\lambda)','Digits of Precision in Translation')
mycbar('L')

storage_Mbytes = (2*L+1).*(L+1).*16/1e6;
h8 = figure(8),clf,contour(d,-log10(prec),storage_Mbytes,20),colorbar
set(gca,'ydir','normal')
myplot('Field Storage for Group Translations: 16-byte-cmplx*(2L+1)*(L+1)','Group Dimension (\lambda)','Digits of Precision in Translation')
mycbar('MB')

%%


delta = 200;
t = 200;
nlevs = 1;

l = 1:nlevs;
d = sqrt(2)*delta*(l.^2);

prec = logspace(-1,-12,20);
[D P] = meshgrid(d,prec);
kd = 2*pi*D;
aph = log10(1./P);
L = ceil(kd + 1.8*aph.^(2/3).*kd.^(1/3));
storage_Mbytes = (2*L+1).*(L+1).*16/1e6;


Nl = (t./(delta*l)).^2;

T = sum(storage_Mbytes.*repmat(Nl,length(prec),1),2)

T/1e3

%%



ltop = 100;
nlevs = 8;
i = 1:nlevs;
li = ltop*(1/2).^(i-1);
di = sqrt(3)*li;
prec = 1e-8;
kd = 2*pi*di;
aph = log10(1./prec);
L = ceil(kd + 1.8*aph.^(2/3).*kd.^(1/3));
L = ~mod(L,2) + L;

storage_Mbytes = (2*L+1).*(L+1).*16/1e6;

nb = 4.^(i-1);

S = storage_Mbytes.*nb
sum(S)


%%









%% FMM setup

% create the number of levels and box dimensions
ltop_lam = 10;
Nlevs = 5;
prec = 1e-16;
[li di kd L] = fmmL(ltop_lam,Nlevs,prec);

% create level property structure
Level = fmmLevel(Nlevs,li,di,L);

% create sturcture and precompute intper/filt polynomials
IntFilt = fmmIntFilt(Nlevs,Level);



%% test the inteprolation

l = 3;

f = rand(Level(l+1).I,Level(l+1).J);
f = ones(Level(l+1).I,Level(l+1).J);

tic
f2 = ssfilt(f,Level(l+1).L,Level(l).L);
toc
tic
f3 = fssfilt(f,Level(l+1).L,Level(l).L);
toc
tic
f4 = fmminterp(f,l,Level,IntFilt);
toc

figure(1),clf,imagesc(real(f)),colorbar%,caxis([-1.5 1.5])
figure(2),clf,imagesc(real(f2)),colorbar%,caxis([-1.5 1.5])
figure(3),clf,imagesc(real(f3-f2)),colorbar%,caxis([-1.5 1.5])
figure(4),clf,imagesc(abs(f4)),colorbar
figure(4),clf,imagesc(real(f3-f4)),colorbar


%% test the scalar filter at L-1

filter = 1;
n = 3;

L = Level(n+1).L;

K = Level(n).L;
T2 = Level(n).Theta;
P2 = Level(n).Phi;
tot2 = K^2 + 2*K + 1;


ylm = sphericalY(K,T2,P2,'mono');
rng(3)
%flm1 = rand(tot2,1) + 1i*rand(tot2,1);
flm1 = zeros(tot2,1);

ind1 = lm2ind(L-1,-(L-1),'scalar');
ind2 = lm2ind(L-1,(L-1),'scalar');
grab = ind1:ind2;
flm1(grab) = rand(length(grab),1) + 1i*rand(length(grab),1);

f = sum(repmat(flm1,[1,numel(T2)]).*ylm,1);
f = reshape(f,size(T2));
figure(1),clf,imagesc(real(f)),colorbar


tic
f2 = ssfilt(f,Level(n+1).L,Level(n).L);
toc
tic
f3 = fssfilt(f,Level(n+1).L,Level(n).L);
toc
tic
f4 = fmmfilterLm1(f,n,Level,IntFilt);
toc


figure(1),clf,imagesc(real(f)),colorbar%,caxis([-1.5 1.5])
figure(2),clf,imagesc(real(f2)),colorbar%,caxis([-1.5 1.5])
figure(3),clf,imagesc(real(f3)),colorbar%,caxis([-1.5 1.5])
figure(4),clf,imagesc(real(f2-f3)),colorbar%,caxis([-1.5 1.5])

figure(4),clf,imagesc(real(f4-f3)),colorbar%,caxis([-1.5 1.5])
figure(4),clf,imagesc(real(f4-f2)),colorbar%,caxis([-1.5 1.5])

figure(1),clf,imagesc(imag(f)),colorbar%,caxis([-1.5 1.5])
figure(2),clf,imagesc(imag(f2)),colorbar%,caxis([-1.5 1.5])
figure(3),clf,imagesc(imag(f3)),colorbar%,caxis([-1.5 1.5])
figure(4),clf,imagesc(imag(f4-f3)),colorbar%,caxis([-1.5 1.5])


%%  test the fast vector filter

% 
n = 3;
K = Level(n+1).L;

L = Level(n).L;
%I = 2*L + 1;
%J = L + 1;
%phi = 2*pi*(0:(I-1))/I;
%[muj, wj] = gqwt(J);
muj = Level(n).mu;
wj = Level(n).w;
T1 = Level(n).Theta;
P1 = Level(n).Phi;
tot1 = L^2 + 2*L;

[Bth Bphi Cth Cphi] = BC(L,T1,P1);

blm1 = zeros(tot1,1);
clm1 = zeros(tot1,1);
% %rng(100)
% blm1 = -0.5 - 0.5*1i + rand(tot1,1) + 1i*rand(tot1,1);
% clm1 = -0.5 - 0.5*1i + rand(tot1,1) + 1i*rand(tot1,1);

% blm1 = ones(tot1,1);
% clm1 = ones(tot1,1);


blm1(2) = 0;
clm1(lm2ind(K+1,1)) = 1;

%blm1(lm2ind(K,4)) = 1;
%clm1(lm2ind(K,4)) = 1;
% blm1(lm2ind(K-1,3)) = 1;
% clm1(lm2ind(K-1,3)) = 1;
% blm1(lm2ind(K+1,3)) = 1;
% clm1(lm2ind(K+1,3)) = 1;
% blm1(lm2ind(K-2,1)) = 1;
% clm1(lm2ind(K-2,1)) = 0;

[Fth Fphi] = BCmult(blm1,clm1,Bth,Bphi,Cth,Cphi);
Fth = reshape(Fth,size(T1));
Fphi = reshape(Fphi,size(T1));
%Fphi = zeros(size(Fphi));



% P = 2*K + 1;
% Q = K + 1;
% phi = 2*pi*(0:(P-1))/P;
% [muk, wk] = gqwt(Q);
muk = Level(n+1).mu;
wk = Level(n+1).w;
%theta = acos(muk);
%[T2 P2] = meshgrid(theta,phi);
T2 = Level(n+1).Theta;
P2 = Level(n+1).Phi;

tot2 = K^2 + 2*K;


%fast scalar filter with correction



n = 3;
OO = ones(size(Fth));
[KX KY KZ] = sph2cart(OO,Level(n).Theta,Level(n).Phi);
dx = 0.5*[1 1 1];
shift = exp(1i*2*pi*(KX*dx(1) + KY*dx(2) + KZ*dx(3)));

shift = 1;

%truncatation method
[blm2 clm2] = vst(shift.*Fth,shift.*Fphi,L,muj,wj);
blm2 = blm2(1:tot2);
clm2 = clm2(1:tot2);
[Fth2,Fphi2] = ivst(blm2,clm2,K,muk);


[Fth3, Fphi3] = vstfilterbasic(shift.*Fth,shift.*Fphi,L,muj,wj,K,muk,T1,T2,P2);

[Fth4, Fphi4] = fmmvecfilter(shift.*Fth,shift.*Fphi,n,Level,IntFilt);

[Fth5, Fphi5] = vstfilterbasic2(shift.*Fth,shift.*Fphi,L,muj,wj,K,muk,T1,T2,P2);


% Fth3 = Level(n+1).sinTheta;
% Fphi3 = zeros(size(Fth3));
% OO = ones(size(Fth3));
% [KX KY KZ] = sph2cart(OO,Level(n+1).Theta,Level(n+1).Phi);
% shift = exp(1i*2*pi*(KX*dx(1) + KY*dx(2) + KZ*dx(3)));
% Fth3 = shift.*Fth3;
% Fphi3 = shift.*Fphi3;


figure(1),clf,imagesc(real(Fth)),colorbar;
figure(2),clf,imagesc(imag(Fth)),colorbar;
figure(3),clf,imagesc(real(Fphi)),colorbar;
figure(4),clf,imagesc(imag(Fphi)),colorbar;

figure(5),clf,imagesc(real(Fth2)),colorbar;
figure(6),clf,imagesc(imag(Fth2)),colorbar;
figure(7),clf,imagesc(real(Fphi2)),colorbar;
figure(8),clf,imagesc(imag(Fphi2)),colorbar;


figure(5),clf,imagesc(real(Fth5)),colorbar;
figure(6),clf,imagesc(imag(Fth5)),colorbar;
figure(7),clf,imagesc(real(Fphi5)),colorbar;
figure(8),clf,imagesc(imag(Fphi5)),colorbar;


figure(9),clf,imagesc(real(Fth3)),colorbar;
figure(10),clf,imagesc(imag(Fth3)),colorbar;
figure(11),clf,imagesc(real(Fphi3)),colorbar;
figure(12),clf,imagesc(imag(Fphi3)),colorbar;


figure(13),clf,imagesc(real(Fth4)),colorbar;
figure(14),clf,imagesc(imag(Fth4)),colorbar;
figure(15),clf,imagesc(real(Fphi4)),colorbar;
figure(16),clf,imagesc(imag(Fphi4)),colorbar;

figure(13),clf,imagesc(real(Fth5-Fth3)),colorbar;
figure(14),clf,imagesc(imag(Fth5-Fth3)),colorbar;
figure(15),clf,imagesc(real(Fphi5-Fphi3)),colorbar;
figure(16),clf,imagesc(imag(Fphi5-Fphi3)),colorbar;

% figure(13),clf,imagesc(real(Fth4-Fth3)),colorbar;
% figure(14),clf,imagesc(imag(Fth4-Fth3)),colorbar;
% figure(15),clf,imagesc(real(Fphi4-Fphi3)),colorbar;
% figure(16),clf,imagesc(imag(Fphi4-Fphi3)),colorbar;
% 

[blm3 clm3] = vst(Fth3,Fphi3,K,muk,wk);
[blm5 clm5] = vst(Fth5,Fphi5,K,muk,wk);

[clm3(1:10) clm5(1:10)]
[blm3(1:10) blm5(1:10)]


%% octree


ltop_lam = 10;

if 1
Npts = 100;
pts = ltop_lam*(rand(Npts,3)-0.5);
pts(:,3) = pts(:,3)*0.1 + 0.2*(pts(:,1) - 1)  + 0.2*(pts(:,2) - 1);

end

if 0
Npts = 300;
pts = ltop_lam*(rand(Npts,3)-0.5);
pts(:,3) = pts(:,3)*0.1 + 0.2*(pts(:,1) - 1)  + 0.2*(pts(:,2) - 1);

pts2 = pts;

pts2(:,3) = pts2(:,3)*0.1 - 0.2*(pts2(:,1) - 1)  - 0.2*(pts2(:,2) - 1) - 10;

pts = [pts; pts2];
Npts = length(pts(:,1));

end


if 0
x = ltop_lam*linspace(-.5,.5,10);
y = x;
[X Y] = meshgrid(x,y);
pts = [X(:) Y(:) ones(size(X(:)))];
Npts = length(X(:));
end


if 0
    
    
ltop_lam = 2;

Npts = 2;
dl = 0.5;
x = dl*[1; 1];
y = d1*[-1; 1];
z  = dl*[1; 1];
pts = [x y z];


offset = 4;
pts2 = pts1;
pts2(:,2) = pts2(:,2) + offset;

end


xmin = 0.5*ltop_lam*[-1 -1 -1];
xmax = 0.5*ltop_lam*[1 1 1];

Tree = BuildOctreeMod(pts(:,1),pts(:,2),pts(:,3),1,xmin,xmax);

%%

h9 = figure(1),clf,hold on

scatter3(pts(:,1),pts(:,2),pts(:,3),'filled')

for n=1:length(Tree),
for m=1:length(Tree(n).group)
    cen = Tree(n).group(m).groupcenter;
    len = Tree(n).group(m).cubelength;
    start = cen - len/2;
    ll = len*[1 1 1];
    plotcube(ll,start,0.01*n,[1 0 0]);
end
end
hold off
axis([xmin(1) xmax(1) xmin(2) xmax(2) xmin(3) xmax(3)])
axis equal
myplot2('Octree','x (\lambda)','y (\lambda)','z (\lambda)');





%% FMM setup

% create the number of levels and box dimensions
Nlevs = length(Tree);
prec = 1e-8;
[li di kd L] = fmmL(ltop_lam,Nlevs+1,prec);

% create level property structure
Level = fmmLevel(Nlevs+1,li,di,L);

% create sturcture and precompute intper/filt polynomials
IntFilt = fmmIntFilt(Nlevs+1,Level);

% add field storage for each group
[Tree] = fmmInitializeTree(Tree,Level);

% create level property structure
Nlevs_direct = 2;
LD = fmmLevel(Nlevs_direct,li([1,end]),di([1,end]),L([1,end]));

% create sturcture and precompute intper/filt polynomials
IntFilt_direct = fmmIntFilt(Nlevs_direct,LD);


%% aggregate direct, interp then shift


Z1 = zeros(LD(1).I,LD(1).J);
O1 = ones(size(Z1));

[KX KY KZ] = sph2cart(O1,LD(1).Theta,LD(1).Phi);

if 0
    figure(2),clf,plot3(KX(:),KY(:),KZ(:),'o');
end


Sth = ones(LD(2).I,LD(2).J,Npts);
Sphi = ones(LD(2).I,LD(2).J,Npts);

[Bth Bphi Cth Cphi] = BC(LD(2).L,LD(2).Theta,LD(2).Phi,'norm');
tot = LD(2).L^2 + 2*LD(2).L;
blm1 = zeros(tot,1);
clm1 = zeros(tot,1);
blm1(2) = 1;
clm1(2) = 1;
[Sth Sphi] = BCmult(blm1,clm1,Bth,Bphi,Cth,Cphi);
Sth = reshape(Sth,size(LD(2).Theta));
Sphi = reshape(Sphi,size(LD(2).Theta));

Sth = repmat(Sth,[1,1,Npts]);
Sphi = repmat(Sphi,[1,1,Npts]);


cen = Tree(1).group(1).groupcenter;

Fth = Z1;
Fphi = Z1;

% amps1 = rand(Npts,1) + 1i*rand(Npts,1);
% amps2 = rand(Npts,1) + 1i*rand(Npts,1);

amps1 = ones(Npts,1);
amps2 = ones(Npts,1);

for n=1:Npts,
    dx = cen - pts(n,:);
    % intepolate, shift, and sum
    shift = exp(1i*2*pi*(KX*dx(1) + KY*dx(2) + KZ*dx(3)));
    Fth = Fth + shift.*fmminterp(amps1(n)*Sth(:,:,n),1,LD,IntFilt_direct);
    Fphi = Fphi + shift.*fmminterp(amps2(n)*Sphi(:,:,n),1,LD,IntFilt_direct);   
end


figure(3),clf,imagesc(LD(1).theta,LD(1).phi,real(Fth)),colorbar;
figure(4),clf,imagesc(LD(1).theta,LD(1).phi,real(Fphi)),colorbar;

figure(3),clf,imagesc(LD(1).theta,LD(1).phi,imag(Fth)),colorbar;
figure(4),clf,imagesc(LD(1).theta,LD(1).phi,imag(Fphi)),colorbar;


%% aggregate direct, shift then interp


Sth = ones(LD(2).I,LD(2).J,Npts);
Sphi = ones(LD(2).I,LD(2).J,Npts);
cen = Tree(1).group(1).groupcenter;


% amps1 = rand(Npts,1) + 1i*rand(Npts,1);
% amps2 = rand(Npts,1) + 1i*rand(Npts,1);

amps1 = ones(Npts,1);
amps2 = ones(Npts,1);

Z1 = zeros(LD(2).I,LD(2).J);
O2 = ones(size(Z1));
Fth = Z1;
Fphi = Z1;

[KX KY KZ] = sph2cart(O2,LD(2).Theta,LD(2).Phi);
for n=1:Npts,
    dx = cen - pts(n,:);
    % intepolate, shift, and sum
    shift = exp(1i*2*pi*(KX*dx(1) + KY*dx(2) + KZ*dx(3)));
    Fth = Fth + shift.*Sth(:,:,n);
    Fphi = Fphi + shift.*Sphi(:,:,n);
end


Fth = fmminterp(Fth,1,LD,IntFilt_direct);
Fphi = fmminterp(Fphi,1,LD,IntFilt_direct);   
    
    
figure(3),clf,imagesc(LD(1).theta,LD(1).phi,real(Fth)),colorbar;
figure(4),clf,imagesc(LD(1).theta,LD(1).phi,real(Fphi)),colorbar;

figure(3),clf,imagesc(LD(1).theta,LD(1).phi,imag(Fth)),colorbar;
figure(4),clf,imagesc(LD(1).theta,LD(1).phi,imag(Fphi)),colorbar;


%% aggregate tree

% initialize
[Tree] = fmmInitializeTree(Tree,Level);

[Bth Bphi Cth Cphi] = BC(LD(2).L,LD(2).Theta,LD(2).Phi,'norm');
tot = LD(2).L^2 + 2*LD(2).L;
blm1 = zeros(tot,1);
clm1 = zeros(tot,1);
blm1(2) = 1;
clm1(2) = 1;
[Sth Sphi] = BCmult(blm1,clm1,Bth,Bphi,Cth,Cphi);
Sth = reshape(Sth,size(LD(2).Theta));
Sphi = reshape(Sphi,size(LD(2).Theta));

Sth = repmat(Sth,[1,1,Npts]);
Sphi = repmat(Sphi,[1,1,Npts]);


Tree = fmmAggregateTree(Tree,Sth,Sphi,pts,Nlevs,Level,IntFilt);




%%

n = 1;
figure(1),clf,imagesc(Level(n).theta,Level(n).phi,real(Fth)),colorbar;
figure(2),clf,imagesc(Level(n).theta,Level(n).phi,real(Fphi)),colorbar;

figure(3),clf,imagesc(Level(n).theta,Level(n).phi,imag(Fth)),colorbar;
figure(4),clf,imagesc(Level(n).theta,Level(n).phi,imag(Fphi)),colorbar;

figure(5),clf,imagesc(Level(n).theta,Level(n).phi,abs(Fth)),colorbar;
figure(6),clf,imagesc(Level(n).theta,Level(n).phi,abs(Fphi)),colorbar;

figure(3),clf,imagesc(Level(n).theta,Level(n).phi,real(Tree(n).group(1).Fth)),colorbar;
figure(4),clf,imagesc(Level(n).theta,Level(n).phi,real(Tree(n).group(1).Fphi)),colorbar;

figure(3),clf,imagesc(Level(n).theta,Level(n).phi,real(Tree(n).group(1).Fth - Fth)),colorbar;
figure(4),clf,imagesc(Level(n).theta,Level(n).phi,real(Tree(n).group(1).Fphi - Fphi)),colorbar;

figure(5),clf,imagesc(Level(n).theta,Level(n).phi,imag(Tree(n).group(1).Fth - Fth)),colorbar;
figure(6),clf,imagesc(Level(n).theta,Level(n).phi,imag(Tree(n).group(1).Fphi - Fphi)),colorbar;


%% 


F1 = Fth;
F2 = Fphi;

F1 = Level(1).sinTheta;
F2 = zeros(size(F1));

F1 = ones(size(Level(1).Theta));
F2 = zeros(size(F1));


Sth = ones(LD(1).I,LD(1).J,Npts);
Sphi = ones(LD(1).I,LD(1).J,Npts);

[Bth Bphi Cth Cphi] = BC(LD(1).L,LD(1).Theta,LD(1).Phi);
tot = LD(1).L^2 + 2*LD(1).L;
blm1 = zeros(tot,1);
clm1 = zeros(tot,1);
blm1(10) = 1;
clm1(10) = 1;
[Sth Sphi] = BCmult(blm1,clm1,Bth,Bphi,Cth,Cphi);
F1 = reshape(Sth,size(LD(1).Theta));
F2 = reshape(Sphi,size(LD(1).Theta));



%% test disaggregate direct, filter then shift

cen = Tree(1).group(1).groupcenter;

F1 = Tree(1).group(1).Fth;
F2 = Tree(1).group(1).Fphi;

[Fthd Fphid] = fmmDisaggregateDirect(F1,F2,pts,cen,LD,IntFilt_direct);


%% test disaggregate direct, shift then filter


Z1 = zeros(LD(1).I,LD(1).J);
O2 = zeros(LD(2).I,LD(2).J,Npts);
Fthd2 = O2;
Fphid2 = O2;
[KX KY KZ] = sph2cart(ones(LD(1).I,LD(1).J),LD(1).Theta,LD(1).Phi);

cen = Tree(1).group(1).groupcenter;
for n=1:Npts,
    
    dx = cen - pts(n,:);
    dx = -dx;
    % intepolate, shift, and sum
    shift = exp(2*pi*1i*(KX*dx(1) + KY*dx(2) + KZ*dx(3)));
    [Sth Sphi] = fmmvecfilter(shift.*F1,shift.*F2,1,LD,IntFilt_direct);

    Fthd2(:,:,n) = Sth;
    Fphid2(:,:,n) = Sphi;
    
end





%%


for n=1:Npts,

    if 0
    figure(1),clf,imagesc(LD(2).theta,LD(2).phi,real(Fthd(:,:,n))),colorbar;
    figure(2),clf,imagesc(LD(2).theta,LD(2).phi,real(Fphid(:,:,n))),colorbar;
 
    figure(3),clf,imagesc(LD(2).theta,LD(2).phi,real(Fthd2(:,:,n))),colorbar;
    figure(4),clf,imagesc(LD(2).theta,LD(2).phi,real(Fphid2(:,:,n))),colorbar;
    end
    
    if 1
    figure(1),clf,imagesc(LD(2).theta,LD(2).phi,abs(Fthd(:,:,n))),colorbar;
    figure(2),clf,imagesc(LD(2).theta,LD(2).phi,abs(Fphid(:,:,n))),colorbar;
 
    figure(3),clf,imagesc(LD(2).theta,LD(2).phi,abs(Fthd2(:,:,n))),colorbar;
    figure(4),clf,imagesc(LD(2).theta,LD(2).phi,abs(Fphid2(:,:,n))),colorbar;
    end
    
    if 0
    
    figure(1),clf,imagesc(LD(2).theta,LD(2).phi,real(Fthd(:,:,n)-Fthd2(:,:,n))),colorbar;
    figure(2),clf,imagesc(LD(2).theta,LD(2).phi,real(Fphid(:,:,n)-Fphid2(:,:,n))),colorbar;
    
    
    figure(3),clf,imagesc(LD(2).theta,LD(2).phi,imag(Fthd(:,:,n)-Fthd2(:,:,n))),colorbar;
    figure(4),clf,imagesc(LD(2).theta,LD(2).phi,imag(Fphid(:,:,n)-Fphid2(:,:,n))),colorbar;
    end
    
    pause
end

%% test disaggregate tree, filter then shift

Tree2 = Tree;

% initialize
for n=1:Nlevs,
   ng = length(Tree2(n).group);
   Z = zeros(Level(n).I,Level(n).J);
   for g=1:ng,
      Tree2(n).group(g).Fth = Z;
      Tree2(n).group(g).Fphi = Z;        
   end
end

% initailize first level
%O1 = ones(Level(1).I,Level(1).J);
Tree2(1).group(1).Fth = F1;
Tree2(1).group(1).Fphi = F2;


% figure(3),clf,imagesc(Level(n).theta,Level(n).phi,real(Tree(n).group(1).Fth)),colorbar;
% figure(4),clf,imagesc(Level(n).theta,Level(n).phi,real(Tree(n).group(1).Fphi)),colorbar;

% filter then shift
% levels Nlevs-1 to 1

for n=1:(Nlevs-1),
for g=1:length(Tree2(n).group),
    
    OO = ones(Level(n+1).I,Level(n+1).J);
    [KX KY KZ] = sph2cart(OO,Level(n+1).Theta,Level(n+1).Phi);
    
    % filter the field from the group to level of child
    [Fthfilt Fphifilt] = fmmvecfilter(Tree2(n).group(g).Fth,Tree2(n).group(g).Fphi,n,Level,IntFilt);
    
    if 0
        if 1
        figure(1),clf,imagesc(abs(Tree2(n).group(g).Fth)),colorbar%,caxis([0 1])
        title(num2str([n g]))
        figure(2),clf,imagesc(real(Tree2(n).group(g).Fth)),colorbar%,caxis([-1 1])
        figure(3),clf,imagesc(imag(Tree2(n).group(g).Fth)),colorbar%,caxis([-1 1])

        figure(4),clf,imagesc(abs(Fthfilt)),colorbar%,caxis([0 1])
        figure(5),clf,imagesc(real(Fthfilt)),colorbar%,caxis([-1 1])
        figure(6),clf,imagesc(imag(Fthfilt)),colorbar%,caxis([-1 1])        
                
        pause
        end
        
        if 0
        figure(1),clf,imagesc(real(Tree2(n).group(g).Fth)),colorbar,caxis([0 1])
        figure(2),clf,imagesc(real(Tree2(n).group(g).Fphi)),colorbar,caxis([0 1])
        figure(3),clf,imagesc(real(Fthfilt)),colorbar,caxis([0 1])
        title(num2str([n g]))
        figure(4),clf,imagesc(real(Fphifilt)),colorbar,caxis([0 1])
        pause
        end
        
        if 0
        figure(1),clf,imagesc(imag(Tree2(n).group(g).Fth)),colorbar,caxis([0 1])
        title(num2str([n g]))
        figure(2),clf,imagesc(imag(Tree2(n).group(g).Fphi)),colorbar,caxis([0 1])
        figure(3),clf,imagesc(imag(Fthfilt)),colorbar,caxis([0 1])
        title(num2str([n g]))
        figure(4),clf,imagesc(imag(Fphifilt)),colorbar,caxis([0 1])
        pause
        end
        
        if 0        
        figure(1),clf,imagesc(angle(Tree2(n).group(g).Fth)),colorbar,caxis([0 1])
        title(num2str([n g]))
        figure(2),clf,imagesc(angle(Tree2(n).group(g).Fphi)),colorbar,caxis([0 1])
        figure(3),clf,imagesc(angle(Fthfilt)),colorbar,caxis([0 1])
        title(num2str([n g]))
        figure(4),clf,imagesc(angle(Fphifilt)),colorbar,caxis([0 1])
        pause
        end
    end
    % shift to child
    for c=1:length(Tree2(n).group(g).child),
        child = Tree2(n).group(g).child(c);
        cen = Tree2(n).group(g).groupcenter;
        xx = Tree2(n+1).group(child).groupcenter;
        dx = cen - xx;
        dx = -dx;
        
        shift = exp(1i*2*pi*(KX*dx(1) + KY*dx(2) + KZ*dx(3)));
        Tree2(n+1).group(child).Fth = shift.*Fthfilt;
        Tree2(n+1).group(child).Fphi = shift.*Fphifilt;

        if 0
            figure(1),clf,hold on

            scatter3(pts(:,1),pts(:,2),pts(:,3),'filled')
            scatter3(cen(1),cen(2),cen(3),'k','filled')
            scatter3(xx(1),xx(2),xx(3),'g','filled')
            quiver3(cen(1),cen(2),cen(3),dx(1),dx(2),dx(3))


            for nn=1:length(Tree),
            for m=1:length(Tree(nn).group)
                cen = Tree(nn).group(m).groupcenter;
                len = Tree(nn).group(m).cubelength;
                start = cen - len/2;
                ll = len*[1 1 1];
                plotcube(ll,start,0.00*nn,[1 0 0]);
            end
            end
            hold off
            axis([xmin(1) xmax(1) xmin(2) xmax(2) xmin(3) xmax(3)])
            axis equal
            pause
        end
    end
            
        
    end

end


Z2 = zeros(Level(Nlevs+1).I,Level(Nlevs+1).J,Npts);
Fth = Z2;
Fphi = Z2;

n=Nlevs;
for g=1:length(Tree2(n).group),
    OO = ones(Level(n+1).I,Level(n+1).J);
    [KX KY KZ] = sph2cart(OO,Level(n+1).Theta,Level(n+1).Phi);
    
    [Fthfilt Fphifilt] = fmmvecfilter(Tree2(n).group(g).Fth,Tree2(n).group(g).Fphi,n,Level,IntFilt);

    for c=1:length(Tree2(n).group(g).child)
        
        child = Tree2(n).group(g).child(c);
        cen = Tree2(n).group(g).groupcenter;
        xx = pts(child,:);
        dx = cen - xx;
        dx = -dx;
        shift = exp(1i*2*pi*(KX*dx(1) + KY*dx(2) + KZ*dx(3)));
        
        Fth(:,:,child) = shift.*Fthfilt;
        Fphi(:,:,child) = shift.*Fphifilt;
    end
end



%% shift then filter

Tree2 = Tree;

% initialize
for n=1:Nlevs,
   ng = length(Tree2(n).group);
   Z = zeros(Level(n).I,Level(n).J);
   for g=1:ng,
      Tree2(n).group(g).Fth = Z;
      Tree2(n).group(g).Fphi = Z;        
   end
end

% initailize first level
%O1 = ones(Level(1).I,Level(1).J);
Tree2(1).group(1).Fth = F1;
Tree2(1).group(1).Fphi = F2;

for n=1:(Nlevs-1),
for g=1:length(Tree2(n).group),
    
    OO = ones(Level(n).I,Level(n).J);
    [KX KY KZ] = sph2cart(OO,Level(n).Theta,Level(n).Phi);
    
    for c=1:length(Tree2(n).group(g).child)
        
        child = Tree(n).group(g).child(c);
        dx = Tree2(n).group(g).groupcenter - Tree2(n+1).group(child).groupcenter;
        dx = -dx;
        shift = exp(1i*2*pi*(KX*dx(1) + KY*dx(2) + KZ*dx(3)));
        [Fthfilt Fphifilt] = fmmvecfilter(shift.*Tree2(n).group(g).Fth,shift.*Tree2(n).group(g).Fphi,n,Level,IntFilt);

        Tree2(n+1).group(child).Fth = Fthfilt;
        Tree2(n+1).group(child).Fphi = Fphifilt;
    end
end
end


Z2 = zeros(Level(Nlevs+1).I,Level(Nlevs+1).J,Npts);
Fth = Z2;
Fphi = Z2;

n=Nlevs;
for g=1:length(Tree2(n).group),
    OO = ones(Level(n).I,Level(n).J);
    [KX KY KZ] = sph2cart(OO,Level(n).Theta,Level(n).Phi);

    for c=1:length(Tree2(n).group(g).child)
        
        child = Tree2(n).group(g).child(c);
        dx = Tree2(n).group(g).groupcenter - pts(child,:);
        dx = -dx;
        shift = exp(1i*2*pi*(KX*dx(1) + KY*dx(2) + KZ*dx(3)));
        [Fthfilt Fphifilt] = fmmvecfilter(shift.*Tree2(n).group(g).Fth,shift.*Tree2(n).group(g).Fphi,n,Level,IntFilt);

        Fth(:,:,child) = Fthfilt;
        Fphi(:,:,child) = Fphifilt;
    end
end

%%

F1 = Tree(1).group(1).Fth;
F2 = Tree(1).group(1).Fphi;

Tree2 = Tree;
[Fth Fphi] = fmmDisaggregateTree1(F1,F2,Tree2,pts,Nlevs,Level,IntFilt);


%%


for n=1:Npts,
    if 0
    figure(1),clf,imagesc(LD(2).theta,LD(2).phi,real(Fthd(:,:,n))),colorbar;
    figure(2),clf,imagesc(LD(2).theta,LD(2).phi,real(Fphid(:,:,n))),colorbar;
    
    figure(3),clf,imagesc(Level(Nlevs).theta,Level(Nlevs).phi,real(Fth(:,:,n))),colorbar;
    title(num2str(norm(pts(n,:))))
    figure(4),clf,imagesc(Level(Nlevs).theta,Level(Nlevs).phi,real(Fphi(:,:,n))),colorbar
    end
    
    if 0
    figure(1),clf,imagesc(LD(2).theta,LD(2).phi,imag(Fthd(:,:,n))),colorbar;
    figure(2),clf,imagesc(LD(2).theta,LD(2).phi,imag(Fphid(:,:,n))),colorbar;
    
    figure(3),clf,imagesc(Level(Nlevs).theta,Level(Nlevs).phi,imag(Fth(:,:,n))),colorbar;
    figure(4),clf,imagesc(Level(Nlevs).theta,Level(Nlevs).phi,imag(Fphi(:,:,n))),colorbar
    end
    
    if 0
    figure(1),clf,imagesc(LD(2).theta,LD(2).phi,abs(Fthd(:,:,n))),colorbar%,caxis([0.955 1.005]);
    figure(2),clf,imagesc(LD(2).theta,LD(2).phi,abs(Fphid(:,:,n))),colorbar%,caxis([0.955 1.005]);
    
    figure(3),clf,imagesc(Level(Nlevs).theta,Level(Nlevs).phi,abs(Fth(:,:,n))),colorbar%,caxis([0.955 1.005]);
    figure(4),clf,imagesc(Level(Nlevs).theta,Level(Nlevs).phi,abs(Fphi(:,:,n))),colorbar%,caxis([0.955 1.005]);
    end
    
    if 0
    figure(3),clf,imagesc(LD(2).theta,LD(2).phi,angle(Fthd(:,:,n))),colorbar;
    figure(4),clf,imagesc(LD(2).theta,LD(2).phi,angle(Fphid(:,:,n))),colorbar;
    
    figure(5),clf,imagesc(Level(Nlevs).theta,Level(Nlevs).phi,angle(Fth(:,:,n))),colorbar;
    figure(6),clf,imagesc(Level(Nlevs).theta,Level(Nlevs).phi,angle(Fphi(:,:,n))),colorbar
    end
    
    if 1
    figure(3),clf,imagesc(LD(2).theta,LD(2).phi,real(Fth(:,:,n)-Fthd(:,:,n))),colorbar;
    figure(4),clf,imagesc(LD(2).theta,LD(2).phi,real(Fphi(:,:,n)-Fphid(:,:,n))),colorbar;
    
    figure(5),clf,imagesc(LD(2).theta,LD(2).phi,imag(Fth(:,:,n)-Fthd(:,:,n))),colorbar;
    figure(6),clf,imagesc(LD(2).theta,LD(2).phi,imag(Fphi(:,:,n)-Fphid(:,:,n))),colorbar;
    end
    
    

    if 0
        figure(10),clf,hold on

        scatter3(pts(:,1),pts(:,2),pts(:,3),'filled')
        %scatter3(cen(1),cen(2),cen(3),'k','filled')
        scatter3(pts(n,1),pts(n,2),pts(n,3),'g','filled')

%         for nn=1:length(Tree),
%         for m=1:length(Tree(nn).group)
%             cen = Tree(nn).group(m).groupcenter;
%             len = Tree(nn).group(m).cubelength;
%             start = cen - len/2;
%             ll = len*[1 1 1];
%             plotcube(ll,start,0.01*nn,[1 0 0]);
%         end
%         end
        hold off
        axis([xmin(1) xmax(1) xmin(2) xmax(2) xmin(3) xmax(3)])
        axis equal
        view([0 0 1])
    end
    
    pause
end



%%

n = 1;
g = 1;
Fthtmp = Tree2(n).group(g).Fth;
Fphitmp = Tree2(n).group(g).Fphi;
Fphitmp = zeros(size(Fthtmp));
Fthtmp = F1;
Fphitmp = F2;

K = Level(n).L;

m = length(Level);
L = Level(m).L;
muj = Level(m).mu;
wj = Level(m).w;
T1 = Level(m).Theta;
P1 = Level(m).Phi;
tot1 = L^2 + 2*L;





muk = Level(n).mu;
wk = Level(n).w;
%theta = acos(muk);
%[T2 P2] = meshgrid(theta,phi);
T2 = Level(n).Theta;
P2 = Level(n).Phi;

tot2 = K^2 + 2*K;


[blm1, clm1] = vst(Fthtmp,Fphitmp,K,muk,wk);

% truncatation method
blm2 = blm1(1:tot2);
clm2 = clm1(1:tot2);
[I J] = size(T1);
[Fth2,Fphi2] = ivst(blm2,clm2,L,muj,I,J);

% fast scalar filter with correction

[Fth3, Fphi3] = vstfilterbasic(Fthtmp,Fphitmp,K,muk,wk,L,muj,T2,T1,P1);

[Fth4, Fphi4] = fmmvecfilter(Fthtmp,Fphitmp,n,LD,IntFilt_direct);


figure(1),clf,imagesc(real(Fthtmp)),colorbar;
figure(2),clf,imagesc(imag(Fthtmp)),colorbar;
figure(3),clf,imagesc(real(Fphitmp)),colorbar;
figure(4),clf,imagesc(imag(Fphitmp)),colorbar;

figure(5),clf,imagesc(real(Fth2)),colorbar;
figure(6),clf,imagesc(imag(Fth2)),colorbar;
figure(7),clf,imagesc(real(Fphi2)),colorbar;
figure(8),clf,imagesc(imag(Fphi2)),colorbar;



figure(1),clf,imagesc(abs(Fthtmp)),colorbar;
figure(3),clf,imagesc(abs(Fphitmp)),colorbar;
figure(5),clf,imagesc(abs(Fth2)),colorbar;
figure(7),clf,imagesc(abs(Fphi2)),colorbar;


figure(9),clf,imagesc(real(Fth3)),colorbar;
figure(10),clf,imagesc(imag(Fth3)),colorbar;
figure(11),clf,imagesc(real(Fphi3)),colorbar;
figure(12),clf,imagesc(imag(Fphi3)),colorbar;

figure(13),clf,imagesc(real(Fth4)),colorbar;
figure(14),clf,imagesc(imag(Fth4)),colorbar;
figure(15),clf,imagesc(real(Fphi4)),colorbar;
figure(16),clf,imagesc(imag(Fphi4)),colorbar;

figure(13),clf,imagesc(real(Fth2-Fth3)),colorbar;
figure(14),clf,imagesc(imag(Fth2-Fth3)),colorbar;
figure(15),clf,imagesc(real(Fphi2-Fphi3)),colorbar;
figure(16),clf,imagesc(imag(Fphi2-Fphi3)),colorbar;

figure(13),clf,imagesc(real(Fth4-Fth3)),colorbar;
figure(14),clf,imagesc(imag(Fth4-Fth3)),colorbar;
figure(15),clf,imagesc(real(Fphi4-Fphi3)),colorbar;
figure(16),clf,imagesc(imag(Fphi4-Fphi3)),colorbar;





%% fast vector filter


% create the number of levels and box dimensions
ltop_lam = 1;
Nlevs = 2;
prec = 1e-8;
[li di kd L] = fmmL(ltop_lam,Nlevs,prec);


% create level property structure
Level = fmmLevel(Nlevs,li,di,L);

% create sturcture and precompute intper/filt polynomials
IntFilt = fmmIntFilt(Nlevs,Level);


%%
% 

n = 1;
K = Level(n).L;
%I = 2*L + 1;
%J = L + 1;
%phi = 2*pi*(0:(I-1))/I;
%[muj, wj] = gqwt(J);
muk = Level(n).mu;
wk = Level(n).w;
T2 = Level(n).Theta;
P2 = Level(n).Phi;
tot2 = K^2 + 2*K;


% spectrum of the shift
O1 = ones(Level(n).I,Level(n).J);
[KX KY KZ] = sph2cart(O1,Level(n).Theta,Level(n).Phi);
x = 1*di(n)/2*[1 1 1];
shift = exp(1i*2*pi*(KX*x(1)+KY*x(2)+KZ*x(3)));
Fth = shift.*ones(size(T2));
Fphi = shift.*zeros(size(T2));


[Bth Bphi Cth Cphi] = BC(K,T2,P2);
[Fth2 Fphi2] = BCmult(blm1,clm1,Bth,Bphi,Cth,Cphi);
Fth2 = reshape(Fth2,size(T2));
Fphi2 = reshape(Fphi2,size(T2));

[blm2, clm2] = vst(Fth2,zeros(size(Fphi2)),K,muk,wk);

figure(1),clf,hold all
plot(abs([blm1 clm1]))
plot(abs([blm2 clm2]))
hold off

Fth3 = Fth;
Fphi3 = Fphi;
for n=1:20,

[blm2, clm2] = vst(Fth3,zeros(size(Fphi3)),K,muk,wk);

    
[Bth Bphi Cth Cphi] = BC(K,T2,P2);
[Fth3 Fphi3] = BCmult(blm2,clm2,Bth,Bphi,Cth,Cphi);
Fth3 = reshape(Fth3,size(T2));
Fphi3 = reshape(Fphi3,size(T2));

end


figure(1),clf,imagesc(real(Fth)),colorbar;
figure(2),clf,imagesc(imag(Fth)),colorbar;
figure(3),clf,imagesc(real(Fphi)),colorbar;
figure(4),clf,imagesc(imag(Fphi)),colorbar;

figure(5),clf,imagesc(real(Fth2)),colorbar;
figure(6),clf,imagesc(imag(Fth2)),colorbar;
figure(7),clf,imagesc(real(Fphi2)),colorbar;
figure(8),clf,imagesc(imag(Fphi2)),colorbar;

figure(9),clf,imagesc(real(Fth3)),colorbar;
figure(10),clf,imagesc(imag(Fth3)),colorbar;
figure(11),clf,imagesc(real(Fphi3)),colorbar;
figure(12),clf,imagesc(imag(Fphi3)),colorbar;

figure(13),clf,imagesc(real(Fth3 - Fth)),colorbar;
figure(14),clf,imagesc(imag(Fth3 - Fth)),colorbar;
figure(15),clf,imagesc(real(Fphi3- Fphi)),colorbar;
figure(16),clf,imagesc(imag(Fphi3- Fphi)),colorbar;






%



[Bth Bphi Cth Cphi] = BC(K,T2,P2);
[Fth Fphi] = BCmult(blm1,clm1,Bth,Bphi,Cth,Cphi);
Fth = reshape(Fth,size(T2));
Fphi = reshape(Fphi,size(T2));



[blm1, clm1] = vst(F1,F2,Level(1).L,Level(1).mu,Level(1).w);


plot(abs([blm1 clm1]))
[Bth Bphi Cth Cphi] = BC(K,T2,P2);
[Fth Fphi] = BCmult(blm1,clm1,Bth,Bphi,Cth,Cphi);
Fth = reshape(Fth,size(T2));
Fphi = reshape(Fphi,size(T2));


figure(1),clf,imagesc(real(Fth)),colorbar;
figure(2),clf,imagesc(imag(Fth)),colorbar;
figure(3),clf,imagesc(real(Fphi)),colorbar;
figure(4),clf,imagesc(imag(Fphi)),colorbar;

[blm2, clm2] = vst(zeros(size(Fth)),Fphi,K,muk,wk);

plot(abs([blm2 clm2]))






figure(5),clf,imagesc(real(F1)),colorbar;
figure(6),clf,imagesc(imag(F1)),colorbar;
figure(7),clf,imagesc(real(F2)),colorbar;
figure(8),clf,imagesc(imag(F2)),colorbar;



n = 2;
L = Level(n).L;
muj = Level(n).mu;
wj = Level(n).w;
T1 = Level(n).Theta;
P1 = Level(n).Phi;
tot1 = L^2 + 2*L;

blm = blm1(1:tot1);
clm = clm1(1:tot1);
[Bth Bphi Cth Cphi] = BC(L,T1,P1);
[Fth2 Fphi2] = BCmult(blm,clm,Bth,Bphi,Cth,Cphi);
Fth2 = reshape(Fth2,size(T1));
Fphi2 = reshape(Fphi2,size(T1));


figure(5),clf,imagesc(real(Fth2)),colorbar;
figure(6),clf,imagesc(imag(Fth2)),colorbar;
figure(7),clf,imagesc(real(Fphi2)),colorbar;
figure(8),clf,imagesc(imag(Fphi2)),colorbar;



[Fth3, Fphi3] = fmmvecfilter(F1,F2,1,Level,IntFilt);


% figure(1),clf,imagesc(real(Fth)),colorbar;
% figure(2),clf,imagesc(imag(Fth)),colorbar;
% figure(3),clf,imagesc(real(Fphi)),colorbar;
% figure(4),clf,imagesc(imag(Fphi)),colorbar;
% 
% figure(5),clf,imagesc(real(Fth2)),colorbar;
% figure(6),clf,imagesc(imag(Fth2)),colorbar;
% figure(7),clf,imagesc(real(Fphi2)),colorbar;
% figure(8),clf,imagesc(imag(Fphi2)),colorbar;



figure(9),clf,imagesc(real(Fth3)),colorbar;
figure(10),clf,imagesc(imag(Fth3)),colorbar;
figure(11),clf,imagesc(real(Fphi3)),colorbar;
figure(12),clf,imagesc(imag(Fphi3)),colorbar;




% figure(1),clf,imagesc(abs(Fth)),colorbar;
% figure(2),clf,imagesc(abs(Fphi)),colorbar;
% figure(3),clf,imagesc(abs(Fth4)),colorbar;
% figure(4),clf,imagesc(abs(Fphi4)),colorbar;

figure(13),clf,imagesc(real(Fth2-Fth3)),colorbar;
figure(14),clf,imagesc(imag(Fth2-Fth3)),colorbar;
figure(15),clf,imagesc(real(Fphi2-Fphi3)),colorbar;
figure(16),clf,imagesc(imag(Fphi2-Fphi3)),colorbar;

% figure(13),clf,imagesc(real(Fth4-Fth3)),colorbar;
% figure(14),clf,imagesc(imag(Fth4-Fth3)),colorbar;
% figure(15),clf,imagesc(real(Fphi4-Fphi3)),colorbar;
% figure(16),clf,imagesc(imag(Fphi4-Fphi3)),colorbar;
% 




%%
% 

n = 1;
K = Level(n).L;
%I = 2*L + 1;
%J = L + 1;
%phi = 2*pi*(0:(I-1))/I;
%[muj, wj] = gqwt(J);
muk = Level(n).mu;
wk = Level(n).w;
T2 = Level(n).Theta;
P2 = Level(n).Phi;
tot2 = K^2 + 2*K;

Fth = rand(size(T2));
Fphi = zeros(size(T2));
[blm1, clm1] = vst(Fth,Fphi,K,muk,wk);


error = zeros(2,tot2,4);
for bc = 1:2,
for ind=1:tot2,
    mydisp(ind,tot2)
    
    
blm1 = zeros(tot2,1);
clm1 = zeros(tot2,1);
if bc
    blm1(ind) = 1;
else
    clm1(ind) = 1;
end

% plot(abs([blm1 clm1]))
[Bth Bphi Cth Cphi] = BC(K,T2,P2);
[Fth Fphi] = BCmult(blm1,clm1,Bth,Bphi,Cth,Cphi);
Fth = reshape(Fth,size(T2));
Fphi = reshape(Fphi,size(T2));


if 1
figure(1),clf,imagesc(real(Fth)),colorbar;
figure(2),clf,imagesc(imag(Fth)),colorbar;
figure(3),clf,imagesc(real(Fphi)),colorbar;
figure(4),clf,imagesc(imag(Fphi)),colorbar;
end


n = 2;
L = Level(n).L;
muj = Level(n).mu;
wj = Level(n).w;
T1 = Level(n).Theta;
P1 = Level(n).Phi;
tot1 = L^2 + 2*L;

blm = blm1(1:tot1);
clm = clm1(1:tot1);
[Bth Bphi Cth Cphi] = BC(L,T1,P1);
[Fth2 Fphi2] = BCmult(blm,clm,Bth,Bphi,Cth,Cphi);
Fth2 = reshape(Fth2,size(T1));
Fphi2 = reshape(Fphi2,size(T1));


if 1
figure(5),clf,imagesc(real(Fth2)),colorbar;
figure(6),clf,imagesc(imag(Fth2)),colorbar;
figure(7),clf,imagesc(real(Fphi2)),colorbar;
figure(8),clf,imagesc(imag(Fphi2)),colorbar;
end


[Fth3, Fphi3] = fmmvecfilter(Fth,Fphi,1,Level,IntFilt);



% figure(1),clf,imagesc(real(Fth)),colorbar;
% figure(2),clf,imagesc(imag(Fth)),colorbar;
% figure(3),clf,imagesc(real(Fphi)),colorbar;
% figure(4),clf,imagesc(imag(Fphi)),colorbar;
% 
% figure(5),clf,imagesc(real(Fth2)),colorbar;
% figure(6),clf,imagesc(imag(Fth2)),colorbar;
% figure(7),clf,imagesc(real(Fphi2)),colorbar;
% figure(8),clf,imagesc(imag(Fphi2)),colorbar;



figure(9),clf,imagesc(real(Fth3)),colorbar;
figure(10),clf,imagesc(imag(Fth3)),colorbar;
figure(11),clf,imagesc(real(Fphi3)),colorbar;
figure(12),clf,imagesc(imag(Fphi3)),colorbar;




% figure(1),clf,imagesc(abs(Fth)),colorbar;
% figure(2),clf,imagesc(abs(Fphi)),colorbar;
% figure(3),clf,imagesc(abs(Fth4)),colorbar;
% figure(4),clf,imagesc(abs(Fphi4)),colorbar;

if 1
figure(13),clf,imagesc(real(Fth2-Fth3)),colorbar;
figure(14),clf,imagesc(imag(Fth2-Fth3)),colorbar;
figure(15),clf,imagesc(real(Fphi2-Fphi3)),colorbar;
figure(16),clf,imagesc(imag(Fphi2-Fphi3)),colorbar;
end

% figure(13),clf,imagesc(real(Fth4-Fth3)),colorbar;
% figure(14),clf,imagesc(imag(Fth4-Fth3)),colorbar;
% figure(15),clf,imagesc(real(Fphi4-Fphi3)),colorbar;
% figure(16),clf,imagesc(imag(Fphi4-Fphi3)),colorbar;
% 

e1 = sum(sum(abs(real(Fth2)-real(Fth3))));
e2 = sum(sum(abs(imag(Fth2)-imag(Fth3))));
e3 = sum(sum(abs(real(Fphi2)-real(Fphi3))));
e4 = sum(sum(abs(imag(Fphi2)-imag(Fphi3))));
tmp = [e1 e2 e3 e4];
error(bc,ind,:) = tmp;

pause
end
end

%%

figure(1),clf,hold all
plot(squeeze(error(1,:,:)))
plot(squeeze(error(2,:,:)))
hold off



