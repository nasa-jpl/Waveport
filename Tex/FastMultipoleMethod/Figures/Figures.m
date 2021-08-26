 
direc = '/Users/mshaynes/Desktop/Book/Waveport/Tex/FastMultipoleMethod/Figures/';

%% Gaussian quadrature

h1 = figure(1),clf,hold all
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

if 0
    saveimage(h1,[direc 'gaussquad'],'epsc');
end


%% Spherical nodes

L = 5;
I = 2*L+1;
J = L+1;

[xj wj] = legpts(J);
theta_j = acos(xj);
phi_i = (0:(I-1))*2*pi/I;

[P T] = meshgrid(phi_i,theta_j);

h1=figure(1),clf
hold on
plot(180*[-1 1],90*[1 1],'k--')
plot([0 0 ],180*[0 1],'k--')
scatter(180/pi*wrap(P(:)),180/pi*T(:),50,'filled')
hold off
axis equal
myplot({'Gaussian Quadratue Sampling on Sphere';['L = ' num2str(L)]},'\phi_i (deg)','\theta_j (deg)')
set(gca,'ydir','reverse')
ylab = [0 45 90 135 180];
yticks(ylab)
yticklabels(num2str(ylab'))
xlab = [-180 -90 0 90 180];
xticks(xlab)
xticklabels(num2str(xlab'))
ylim(180/pi*[0 pi])
xlim(180/pi*pi*[-1 1])


if 0
    saveimage(h1,[direc 'quadsphere' num2str(L)],'epsc');
end




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

h1 = figure(1),clf
imagesc(180/pi*phi,180/pi*th,real(tl))
colorbar
myplot({'FMM Translation Operator, TL ';['L = ' num2str(L), ', kX = ' num2str(kX)]},'\phi (deg)','\theta (deg)')
set(gca,'ydir','reverse')
grid off
ylim(180/pi*[0 pi])
xlim(180/pi*pi*[-1 1])


if 0
    saveimage(h1,[direc 'TLtheta'],'epsc');
end



%% TL interpolator, interpTLprep

L = 4;
r = 10;
k = 2*pi;



% interpolation
s = 5;
p = 3;
M = s*L;
dth = 2*pi/(2*M+1);
thetao = p*dth;
[tlsamp, dth, thetao, M, theta_samp] = interpTLprep(L,k,r,s,p);
theta = linspace(0,pi,8*M);
tlinterp = interpTL(theta,L,p,tlsamp,dth,thetao,M);

% direct samling
tlsamp = TLth(L,k*r,cos(theta_samp));

h1 = figure(1);,clf,hold all
plot(theta,real(tlinterp),'linewidth',2)
plot(theta,imag(tlinterp),'linewidth',2)
scatter(theta_samp,real(tlsamp),'o','filled')
scatter(theta_samp,imag(tlsamp),'o','filled')
xlimits = [-0.5 3.75];
ylimits = max(abs(tlsamp))*1.2*[-1 1];
h = plot(xlimits,[0 0],'k','linewidth',2);
plot([0 0],ylimits,'k--','linewidth',2)
plot([pi pi],ylimits,'k--','linewidth',2)
myplot([],'\theta (radians)','Re/Im T_L(\theta)')
legend('Re T_L(\theta) - Interp','Im T_L(\theta) - Interp',...
    'Re T_L(m\Delta\theta) - Sampled','Im T_L(m\Delta\theta) - Sampled')
xlim(xlimits)
ylim(ylimits)
uistack(h,'down',4)


if 0
    saveimage(h1,[direc 'TLthetainterp'],'epsc');
end


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
    totcutoff = Lmax^2 + 2*Lmax + 1;
    flm((totcutoff+1):end) = 0;
    f = reshape(ylm*flm,size(T2));
end

% function determines whether to filter or interpolate based on size(f),L,K
f2 = ssfilt(f,L,muj,wj,K,muk,wk);

if ~filter_or_interp
    
    h1 = figure(1);clf,imagesc(180/pi*phi,180/pi*theta,real(f)'),colorbar,caxis([-1 3])
    myplot2(['L = ' num2str(L)],'\phi (deg)','\theta (deg)')
    
    h2 = figure(2);clf,imagesc(180/pi*phi2,180/pi*theta2,real(f2)'),colorbar,caxis([-1 3])
    myplot2(['K = ' num2str(K)],'\phi (deg)','\theta (deg)')
    
    if 1
        saveimage(h1,[direc 'filt1'],'epsc');
        saveimage(h2,[direc 'filt2'],'epsc');
    end


else
    
    h3 = figure(3);clf,imagesc(180/pi*phi,180/pi*theta,real(f)'),colorbar,caxis([-2 2])
    myplot2(['K = ' num2str(K)],'\phi (deg)','\theta (deg)')

    h4 = figure(4);clf,imagesc(180/pi*phi2,180/pi*theta2,real(f2)'),colorbar,caxis([-2 2])
    myplot2(['L = ' num2str(L)],'\phi (deg)','\theta (deg)')

    if 1
        saveimage(h3,[direc 'filt3'],'epsc');
        saveimage(h4,[direc 'filt4'],'epsc');
    end
end




%%












%% Old



%% legendreP, legendrePp

Np = 500;

L = 5;
xx = linspace(-1,1,Np);
P = legendreP(L,xx);
plot(P')
Pp = legendrePp(L,xx);

leg = zeros(L+1,Np);
for l = 0:L,
    tmp = legendre(l,xx);
    leg(l+1,:) = tmp(1,:);
end

plot(leg')


f2 = 1/2*(3*xx.^2 - 1);
t1 = P(3,:);
t2 = leg(3,:);

plot([t1'-f2' t2'-f2'])

figure(1),clf,plot([P'])
figure(2),clf,plot(Pp')
Pp2 = diff(P,[],2)/(xx(2)-xx(1));
Pp2 = [Pp2 Pp2(:,end)];
figure(3),clf,plot(Pp2')


%%



%% d/dx tilde Plm

L = 4;
x = linspace(-1,1,100);
plm = Plm(L,x);
dx = x(2)-x(1);

plmp = zeros(size(plm));
for n=1:length(plm(:,1)),
    plmp(n,:) = gradient(plm(n,:))/dx;
    
end


dplm = zeros(size(plm));

for l = 0:L,
for m = -l:l,
    ind = lm2ind(l,m,'scalar');
    if l == 0
        c2 = 0;
    elseif abs(m) <= l-1,        
        ind2 = lm2ind(l-1,m,'scalar');
        c2 = -sqrt(l+1/2)/sqrt(l-1/2)*sqrt((l-m)*(l+m));
        c2 = c2.*plm(ind2,:);
    else
        c2 = 0;
    end
    c1 = l*x.*plm(ind,:);
    dplm(ind,:) = (1./(x.^2-1)).*(c1 + c2);
end
end


figure(1),clf,plot(x,plm')

figure(4),clf,hold all
plot(x,plmp')
plot(x,dplm','--')
hold off



%% d/dtheta tilde Plm(cos(theta))

L = 4;
x = linspace(0,pi,100);
plm = Plm(L,cos(x));
dx = x(2)-x(1);

plmp = zeros(size(plm));
for n=1:length(plm(:,1)),
    plmp(n,:) = gradient(plm(n,:))/dx;
end

dplm = Plmp2(L,x);



figure(1),clf,plot(x,plm')

figure(4),clf,hold all
plot(x,plmp')
plot(x,dplm','--')
hold off


%% d/dtheta tilde Plm(cos(theta))

L = 4;
x = linspace(0,pi,100);
plm = Plm(L,cos(x));
tot = L^2 + 2*L + 1;
[l m] = ind2lm((1:tot),'scalar');
[M S] = ndgrid(m,1./sin(x));
plm = M.*S.*plm;

plm2 = mPlmsin(L,x);

figure(4),clf,hold all
plot(x,plm')
plot(x,plm2','--')
hold off


%%

LL = 60;
running = zeros(LL+1,1);

k = 1;
rji = 400*[ 1 1 1];

for L = 0:LL,

    tl = zeros(size(T));

   % for n=10:10,
   n = 100;
        khat = [kx(n) ky(n) kz(n)];
        tl = TL(L,k,khat,rji);
   % end

    running(L+1) = tl(1);

end

plot(10*log10(abs(running)))


%% TL interpolator

L = 200;
s = 5;
p = 3;
M = s*L;

dth = 2*pi/(2*M+1);
thetao = p*dth;
N = M - L;

r = 500;
k = 2*pi;


[tlsamp, dth, thetao, M, theta_samp] = interpTLprep(L,k,r,s,p);

theta = linspace(0,pi,8*M);
tlinterp = interpTL(theta,L,p,tlsamp,dth,thetao,M);


%
xlimits = [-0.5 3.75];
ylimits = 100*[-0.1 .35];
tlsamp = TLth(L,k,r,theta_samp);
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


%%  Scalar spectral filter

% points on circumference
M = 7;
dth = 2*pi/(2*M+1);
theta = 0:dth:(2*pi-dth);

xt = cos(theta);
yt = sin(theta);
figure(1),clf,plot(xt,yt,'o')
myplot([],[],[])
axis equal

% points on half circumference avoiding poles
N = 7;
dth = pi/(N+1);
theta = dth*((1:(N+1))-1/2) - pi/2;

xt = cos(theta);
yt = sin(theta);
figure(1),clf,plot(xt,yt,'o')
myplot([],[],[])
axis equal


%% standard filter

N = 15;
K = 17;

J = N+1;
I = 2*N+1;
P = K+1;
Q = 2*K+1;

dth = pi/J;
dthprime = pi/P;
dphi = 2*pi/I;
dphiprime = 2*pi/Q;

theta = dth*((1:J)-1/2);
thetaprime = dthprime*((1:P)-1/2) - pi/2;
phi = 0:dphi:(2*pi-dphi);
phiprime = 0:dphiprime:(2*pi-dphiprime);

[P T] = meshgrid(phi,theta);
[Pp Tp] = meshgrid(phiprime,thetaprime);
Npts = numel(T);

% analysis 
tot = N^2 + 2*N + 1;
totprime = K^2 + 2*K + 1;

fnm = zeros(tot,1);
fnm(4) = 1;

ylm = sphericalY(N,T,P,'mono');
f = sum(repmat(fnm,1,Npts).*ylm,1);
f = reshape(f,size(P));

figure(1),clf,imagesc(phi,theta,real(f))
figure(2),clf,imagesc(phi,theta,imag(f))

%
f_theta = sqrt(2*pi)/I*fft(f,[],2);

figure(1),clf,imagesc(phi,theta,abs(f_theta)),colorbar

muj = cos(theta);


%%



lpp = legendrePp(J,muj);
lpp = lpp(end,:);

gqwt = 2./(1-muj.^2)./lpp.^2;

f = muj.^2 + 1;

sum(f.*gqwt)


[xj wt] = lgwt(J,-1,1);

f = xj.^2 + 1;

sum(f.*wt)



%% fmm 1D basic

N = 100;
rng(1);
a = -1;
b = 1;

xj = linspace(a,b,N)'; %2*(rand(N,1)-0.5);
ak = ones(size(xj)); %2*(rand(N,1)-0.5);
xk = xj;

% xj = (2/(b-a))*(xj+a) + 1;
% xk = (2/(b-a))*(xk+a) + 1;
% 
% scale = 2/(b-a);


%xk = sort(xk);

f = zeros(N,1);
ffar = f;
fnear = f;
ffardir = f;

prec = 1e-5;
p = ceil(-logb(5,prec));
s = 2*p;
nlevs = ceil(log2(N/s));

tot = 2^(nlevs+1) - 2;

nbox_nlevs = 2^(nlevs);
width_nlevs = (1/2)^(nlevs-1);
r_nlevs = width_nlevs/2;
cen_nlevs = (0:(nbox_nlevs-1))*width_nlevs - 1 + r_nlevs;

% index the points at the lowest level

xj_ind_nlevs = ceil((xj + 1)/2*nbox_nlevs);
xj_ind_nlevs(xj_ind_nlevs == 0) = 1;
xk_ind_nlevs = ceil((xk + 1)/2*nbox_nlevs);
xk_ind_nlevs(xk_ind_nlevs == 0) = 1;

Phi = zeros(p,tot);
Psi = zeros(p,tot);

% Chebyshev coefficients
ti = cos((2*(1:p)-1)*pi/p/2);

% Compute M, S, T
ML = fmmuj(3*ti./(6+ti),ti);
MR = fmmuj(3*ti./(6-ti),ti);
SL = fmmuj((ti-1)/2,ti);
SR = fmmuj((ti+1)/2,ti);
T1 = fmmuj(3./(ti+6),ti);
T2 = fmmuj(3./(ti+4),ti);
T3 = fmmuj(3./(ti-4),ti);
T4 = fmmuj(3./(ti-6),ti); 


% compute far-field at nlevs boxes, k points
l = nlevs;
for i = 1:(2^l),
    ind = find(xk_ind_nlevs == i);
    for pp=1:p,        
        tmp = ti(pp)./(3*r_nlevs - ti(pp)*(xk(ind)-cen_nlevs(i)));
        Phi(pp,box2ind(l,i)) = sum(ak(ind).*tmp);
    end    
end
 


% far-field at each subinterval, each level
for l=(nlevs-1):-1:2,
   for i=1:(2^l),
       Phi(:,box2ind(l,i)) = ML*Phi(:,box2ind(l+1,2*i-1)) + MR*Phi(:,box2ind(l+1,2*i));
   end    
end

% plot ff at all levels
if 0
    y = linspace(-1,1,N)';
    for l = (nlevs):-1:2,
        B = 2^(l);
        W = (1/2)^(l-1);
        r = W/2;
        cen = (0:(B-1))*W - 1 + r;        
        
        for i = 1:(2^l),
            
            % aggregation
            uu = fmmuj(3*r./(y-cen(i)),ti);
            ffar = uu*Phi(:,box2ind(l,i));
            
            % source points
            xk_ind = ceil((xk + 1)/2*B);
            xk_ind((xk_ind == 0)) = 1;
            indk = find(xk_ind == i);
            
            % ffar exp dir
            phi2 = zeros(p,1);
            for pp=1:p,        
                tmp = ti(pp)./(3*r - ti(pp)*(xk(indk)-cen(i)));
                phi2(pp) = sum(ak(indk).*tmp);
            end    
            
            ffar2 = uu*phi2; 
            
            phi2
            Phi(:,box2ind(l,i))
            
            % direct           
            ffardir = zeros(size(y));
            for k=1:length(indk),
                ffardir = ffardir + (1./(y-xk(indk(k))))*ak(indk(k));
            end
            plot(y,[ffardir ffar]),ylim(1000*[-1 1])
            pause
        end

    end
end


% local expansions each subinterval, each level
for l=1:(nlevs-1),
    for i=1:(2^l),
        if (2*i-3) >= 1
            tmpL1 = T2*Phi(:,box2ind(l+1,2*i-3));
            tmpR1 = T1*Phi(:,box2ind(l+1,2*i-3));
        else
            tmpL1 = 0;
            tmpR1 = 0;
        end
        if (2*i-2) >= 1
            tmp2 = T2*Phi(:,box2ind(l+1,2*i-2));
        else
            tmp2 = 0;
        end
        if (2*i+1) <= 2^(l+1)
            tmp3 = T3*Phi(:,box2ind(l+1,2*i+1));
        else
            tmp3 = 0;
        end
        if (2*i+2) <= 2^(l+1)
            tmpL4 = T4*Phi(:,box2ind(l+1,2*i+2));
            tmpR4 = T3*Phi(:,box2ind(l+1,2*i+2));
        else
            tmpL4 = 0;
            tmpR4 = 0;
        end
        Psi(:,box2ind(l+1,2*i-1)) = SL*Psi(:,box2ind(l,i)) + tmpL1 + tmp3 + tmpL4;
        Psi(:,box2ind(l+1,2*i))   = SR*Psi(:,box2ind(l,i)) + tmpR1 + tmp2 + tmpR4;     
     end
end


l = nlevs;
for i = 1:(2^l),
    % obsevation points
    ind = find(xj_ind_nlevs == i);
    
    % evaluate local expansions at subset of yj
    uu = fmmuj((xj(ind)-cen_nlevs(i))/r_nlevs,ti);
    
    ffar(ind) = uu*Psi(:,box2ind(l,i));

    % compute near terms in neighbor boxes
    ind1 = find(xk_ind_nlevs == i-1);
    ind2 = find(xk_ind_nlevs == i);
    ind3 = find(xk_ind_nlevs == i+1);
    indk = [ind1; ind2; ind3];
    
    for n=1:length(ind),
    for k=1:length(indk),
        if ind(n)~=indk(k)
            fnear(ind(n)) = fnear(ind(n)) + (1./(xj(ind(n))-xk(indk(k))))*ak(indk(k));
        end
    end
    end
    
    indk2 = setxor(indk,1:N);
    for n=1:length(ind),
    for k=1:length(indk2),
        if ind(n)~=indk2(k)
            ffardir(ind(n)) = ffardir(ind(n)) + (1./(xj(ind(n))-xk(indk2(k))))*ak(indk2(k));
        end
    end
    end
    
end

f = fnear + ffar;


%%

N = 1000;
rng(1);
a = -1;
b = 1;
step = 0.01;

xj = linspace(a,b-step,N)'; %2*(rand(N,1)-0.5);
ak = ones(size(xj)); %2*(rand(N,1)-0.5);
xk = xj + step;


% tic 
% [ML,MR,SL,SR,T1,T2,T3,T4,...
%     ti,p,nlevs,r,cen,xj_ind,xk_ind,Phi,Psi,B] = fmm1prep_orig(xk,xj,1e-5);
% toc 
% 
% tic
% [f] = fmm1_orig(ak,xk,xj,ML,MR,SL,SR,T1,T2,T3,T4,...
%             ti,p,nlevs,r,cen,xj_ind,xk_ind,Phi,Psi);
% toc


tic 
[S] = fmm1prep(xk,xj,1e-5);
toc 

tic
[f] = fmm1(ak,S);
toc


% direct
[Xj Xk] = ndgrid(xj,xk);
M = 1./(Xj-Xk);
ind = find(M == Inf);
M(ind) =0;

tic
fdir = M*ak;
toc


%figure(2),clf,plot(xj,[ffardir ffar])
%figure(2),clf,plot(xj,[ffardir ffar])

figure(3),clf,plot(xj,real([fdir-f]))
figure(4),clf,plot(xj,real([fdir f]))

%%




%% fmm 1D SVD

N = 2000;
rng(1);
a = -1;
b = 1;

xj = linspace(a,b,N)'; %2*(rand(N,1)-0.5);
%ak = rand(size(xj)) + 1i*rand(size(xj)); %2*(rand(N,1)-0.5);
ak = ones(size(xj));
xk = xj;

% xj = (2/(b-a))*(xj+a) + 1;
% xk = (2/(b-a))*(xk+a) + 1;
% 
% scale = 2/(b-a);


f = zeros(N,1);
ffar = f;
fnear = f;
ffardir = f;

prec = 1e-6;
p = ceil(-logb(5,prec));
s = 2*p;
nlevs = ceil(log2(N/s));

tot = 2^(nlevs+1) - 2;

nbox_nlevs = 2^(nlevs);
width_nlevs = (1/2)^(nlevs-1);
r_nlevs = width_nlevs/2;
cen_nlevs = (0:(nbox_nlevs-1))*width_nlevs - 1 + r_nlevs;

% index the points at the lowest level

xj_ind_nlevs = ceil((xj + 1)/2*nbox_nlevs);
xj_ind_nlevs(xj_ind_nlevs == 0) = 1;
xk_ind_nlevs = ceil((xk + 1)/2*nbox_nlevs);
xk_ind_nlevs(xk_ind_nlevs == 0) = 1;

% Chebyshev coefficients
ti = cos((2*(1:p)-1)*pi/p/2);
a = 0;
b = 1;
tiL = 0.5*(a+b) + 0.5*(b-a)*cos((2*(1:p)-1)*pi/p/2);
a = -1;
b = 0;
tiR = 0.5*(a+b) + 0.5*(b-a)*cos((2*(1:p)-1)*pi/p/2);

[TjL Ti] = meshgrid(tiL,ti);
[TjR Ti] = meshgrid(tiR,ti);

MML = 1./(3./TjL - Ti);
MMR = 1./(3./TjR - Ti);

phat = rank(MML);
phat= p-1;

[UL EL VL] = svd(MML);
[UR ER VR] = svd(MMR);

UL = UL(:,1:phat);
VL = VL(:,1:phat);
UR = UR(:,1:phat);
VR = VR(:,1:phat);

% Compute M, S, T
ML = fmm1u(3*ti./(6+ti),ti);
MR = fmm1u(3*ti./(6-ti),ti);
SL = fmm1u((ti-1)/2,ti);
SR = fmm1u((ti+1)/2,ti);
T1 = fmm1u(3./(ti+6),ti);
T2 = fmm1u(3./(ti+4),ti);
T3 = fmm1u(3./(ti-4),ti);
T4 = fmm1u(3./(ti-6),ti); 

ML_L = UL.'*ML*UL;
MR_L = UL.'*MR*UL;
ML_R = UR.'*ML*UR;
MR_R = UR.'*MR*UR;
SL_L = VL.'*SL*VL;
SR_L = VL.'*SR*VL;
SL_R = VR.'*SL*VR;
SR_R = VR.'*SR*VR;

TT1 = VL.'*T1*UL;
TT2 = VL.'*T2*UL;
TT3 = VR.'*T3*UR;
TT4 = VR.'*T4*UR;

Phi = zeros(p,tot);
Psi = zeros(p,tot);
PhiL = zeros(phat,tot);
PhiR = zeros(phat,tot);
PsiL = zeros(phat,tot);
PsiR = zeros(phat,tot);



% compute far-field at nlevs boxes, k points
l = nlevs;
for i = 1:(2^l),
    ind = find(xk_ind_nlevs == i);
    for pp=1:p,        
        tmp = ti(pp)./(3*r_nlevs - ti(pp)*(xk(ind)-cen_nlevs(i)));
        Phi(pp,box2ind(l,i)) = sum(ak(ind).*tmp);
    end    
end
 
% compress
l = nlevs;
for i = 1:(2^l),
    ind = box2ind(l,i);
    PhiL(:,ind) = UL.'*Phi(:,ind);
    PhiR(:,ind) = UR.'*Phi(:,ind);
end

% far-field at each subinterval, each level
for l=(nlevs-1):-1:2,
   for i=1:(2^l),
       PhiL(:,box2ind(l,i)) = ML_L*PhiL(:,box2ind(l+1,2*i-1)) + MR_L*PhiL(:,box2ind(l+1,2*i));
       PhiR(:,box2ind(l,i)) = ML_R*PhiR(:,box2ind(l+1,2*i-1)) + MR_R*PhiR(:,box2ind(l+1,2*i));
   end    
end


% local expansions each subinterval, each level
for l=1:(nlevs-1),
    for i=1:(2^l),
        if (2*i-3) >= 1
            c1 = TT2*PhiL(:,box2ind(l+1,2*i-3));
            c2 = TT1*PhiL(:,box2ind(l+1,2*i-3));
        else
            c1 = 0;
            c2 = 0;
        end
        if (2*i-2) >= 1
            c3 = TT2*PhiL(:,box2ind(l+1,2*i-2));
        else
            c3 = 0;
        end
        if (2*i+1) <= 2^(l+1)
            c4 = TT3*PhiR(:,box2ind(l+1,2*i+1));
        else
            c4 = 0;
        end
        if (2*i+2) <= 2^(l+1)
            c5 = TT4*PhiR(:,box2ind(l+1,2*i+2));
            c6 = TT3*PhiR(:,box2ind(l+1,2*i+2));
        else
            c5 = 0;
            c6 = 0;
        end
       % Psi(:,box2ind(l+1,2*i-1)) = SL*Psi(:,box2ind(l,i)) + tmpL1 + tmp3 + tmpL4;
       % Psi(:,box2ind(l+1,2*i))   = SR*Psi(:,box2ind(l,i)) + tmpR1 + tmp2 + tmpR4;
        PsiR(:,box2ind(l+1,2*i-1)) = SL_R*PsiR(:,box2ind(l,i)) + c4 + c5;
        PsiR(:,box2ind(l+1,2*i))   = SR_R*PsiR(:,box2ind(l,i)) + c6;

        PsiL(:,box2ind(l+1,2*i-1)) = SL_L*PsiL(:,box2ind(l,i)) + c1;
        PsiL(:,box2ind(l+1,2*i))   = SR_L*PsiL(:,box2ind(l,i)) + c3 + c2;

     end
end


l = nlevs;
for i = 1:(2^l),
    % obsevation points
    ind = find(xj_ind_nlevs == i);
    
    % evaluate local expansions at subset of yj
    uu = fmm1u((xj(ind)-cen_nlevs(i))/r_nlevs,ti);
    ffar(ind) = uu*(VL*PsiL(:,box2ind(l,i)) + VR*PsiR(:,box2ind(l,i)));

%     uuL = fmm1u((xj(ind)-cen_nlevs(i))/r_nlevs,tiL);
%     uuR = fmm1u((xj(ind)-cen_nlevs(i))/r_nlevs,tiR);
%     ffar(ind) = uuL*VL*PsiL(:,box2ind(l,i)) + uuR*VR*PsiR(:,box2ind(l,i));

    
    % compute near terms in neighbor boxes
    ind1 = find(xk_ind_nlevs == i-1);
    ind2 = find(xk_ind_nlevs == i);
    ind3 = find(xk_ind_nlevs == i+1);
    indk = [ind1; ind2; ind3];
    
    for n=1:length(ind),
        indsum = (xj(ind(n)) ~= xk(indk));
        fnear(ind(n)) = fnear(ind(n)) + sum((1./(xj(ind(n))-xk(indk(indsum)))).*ak(indk(indsum)));
    end    
    
%     indk2 = setxor(indk,1:N);
%     for n=1:length(ind),
%     for k=1:length(indk2),
%         if ind(n)~=indk2(k)
%             ffardir(ind(n)) = ffardir(ind(n)) + (1./(xj(ind(n))-xk(indk2(k))))*ak(indk2(k));
%         end
%     end
%     end
    
end

f = fnear + ffar;




% direct
[Xj Xk] = ndgrid(xj,xk);
M = 1./(Xj-Xk);
ind = find(M == Inf);
M(ind) =0;

fdir = M*ak;



fdir2 = fnear + ffardir;
figure(1),clf,plot(xj,real([fdir fdir2 fnear ffardir]))
figure(2),clf,plot(xj,imag([fdir fdir2 fnear ffardir]))

figure(2),clf,plot(xj,[ffardir ffar])
figure(2),clf,plot(xj,[ffardir ffar])

figure(3),clf,plot(xj,real([fdir f]))
figure(4),clf,plot(xj,imag([fdir f]))

%%

L = 6;
J = L + 1;
I = 2*L + 1;
phi = 2*pi*(0:(I-1))/I;
[muj, wj] = gqwt(J);
theta = acos(muj);

[T P] = meshgrid(theta,phi);
ylm = sphericalY(L,T,P,'mono');

tot = L^2 + 2*L + 1;
rng(1)
flm_o = rand(tot,1) + 1i*rand(tot,1);

f = sum(repmat(flm_o,[1,numel(T)]).*ylm,1);
f = reshape(f,size(T));

flm = sst(f,L,muj,wj);


[l m] = ind2lm((1:tot)','scalar');
[l m flm_o flm]



ff = isst(flm,L,muj);

figure(1),clf,imagesc(real(f)),colorbar
figure(2),clf,imagesc(real(ff)),colorbar


%%

L = 6;
Lsamp = L+1;
J = Lsamp + 1;
I = 2*Lsamp + 1;
phi = 2*pi*(0:(I-1))/I;
[muj, wj] = gqwt(J);
theta = acos(muj);

[T P] = meshgrid(theta,phi);
ylm = sphericalY(L,T,P,'mono');

tot = L^2 + 2*L + 1;
rng(1)
%flm_o = rand(tot,1) + 1i*rand(tot,1);

for ll=0:L,
for m=-ll:ll,
    flm_o = zeros(tot,1);
    flm_o(lm2ind(ll,m,'scalar')) = 1;

    f = sum(repmat(flm_o,[1,numel(T)]).*ylm,1);
    f = reshape(f,size(T));

    flm = sst(f,L,muj,wj);


    [l m] = ind2lm((1:tot)','scalar');
    [l m flm_o flm]

pause
end
end

ff = isst(flm,L,muj);

figure(1),clf,imagesc(real(f)),colorbar
figure(2),clf,imagesc(real(ff)),colorbar

%% test L+1 sampling

L = 5;
Lp = 12;
J = Lp + 1;
I = 2*Lp + 1;
phi = 2*pi*(0:(I-1))/I;
[muj, wj] = gqwt(J);
theta = acos(muj);

[T P] = meshgrid(theta,phi);
ylm = sphericalY(L,T,P,'mono');

tot = L^2 + 2*L + 1;
rng(1)
flm_o = rand(tot,1) + 1i*rand(tot,1);

f = sum(repmat(flm_o,[1,numel(T)]).*ylm,1);
f = reshape(f,size(T));

flm = sst(f,L,muj,wj);


[l m] = ind2lm((1:tot)','scalar');
[l m flm_o flm]



ff = isst(flm,L,muj);

figure(1),clf,imagesc(real(f)),colorbar
figure(2),clf,imagesc(real(ff)),colorbar

%%




flm = sst(f,L,muj,wj);


[l m] = ind2lm((1:tot)','scalar');
[l m flm_o flm]



ff = isst(flm,L,muj);

figure(1),clf,imagesc(real(f)),colorbar
figure(2),clf,imagesc(real(ff)),colorbar





%%
% filter

L = 4;
I = 2*L + 1;
J = L + 1;
phi = 2*pi*(0:(I-1))/I;
[muj, wj] = gqwt(J);
theta = acos(muj);
[T1 P1] = meshgrid(theta,phi);
tot1 = L^2 + 2*L + 1;

ylm = sphericalY(L,T1,P1,'mono');

rng(8)
flm1 = rand(tot1,1) + 1i*rand(tot1,1);
%flm1 = zeros(tot1,1);
%flm1(14) = 1;

f = sum(repmat(flm1,[1,numel(T1)]).*ylm,1);
f = reshape(f,size(T1));
figure(1),clf,imagesc(real(f)),colorbar%,caxis([-1.5 1.5])



K = 20;
P = 2*K + 1;
Q = K + 1;
phi2 = 2*pi*(0:(P-1))/P;
[muk wk] = gqwt(Q);
theta2 = acos(muk);
[T2 P2] = meshgrid(theta2,phi2);
tot2 = K^2 + 2*K + 1;


filt = 0;
if filt
    ylm = sphericalY(K,T2,P2,'mono');
    rng(3)
    %flm1 = rand(tot2,1) + 1i*rand(tot2,1);
    flm1 = zeros(tot2,1);
    
    grab = 6:20;
    flm1(grab) = rand(length(grab),1) + 1i*rand(length(grab),1);

    f = sum(repmat(flm1,[1,numel(T2)]).*ylm,1);
    f = reshape(f,size(T2));
    figure(1),clf,imagesc(real(f)),colorbar
end

% forward

f2 = ssfilt(f,L,K);
[f3] = fssbasic(f,L,K);


figure(2),clf,imagesc(real(f2)),colorbar%,caxis([-1.5 1.5])
figure(3),clf,imagesc(real(f3)),colorbar%,caxis([-1.5 1.5])
figure(4),clf,imagesc(real(f2 - f3)),colorbar

if ~filt
    
    h2 = figure(1);clf,imagesc(180/pi*phi,180/pi*theta,real(f)'),colorbar,caxis([-1 3])
    myplot2(['L = ' num2str(L)],'\phi (deg)','\theta (deg)')
    h2 = myaa(12);
    
    h3 = figure(2);clf,imagesc(180/pi*phi2,180/pi*theta2,real(f2)'),colorbar,caxis([-1 3])
    myplot2(['K = ' num2str(K)],'\phi (deg)','\theta (deg)')
    h3 = myaa(12);

else
    
    h4 = figure(3);clf,imagesc(180/pi*phi,180/pi*theta,real(f)'),colorbar,caxis([-1 2])
    myplot2(['K = ' num2str(K)],'\phi (deg)','\theta (deg)')
        h4 = myaa(12);

    h5 = figure(4);clf,imagesc(180/pi*phi2,180/pi*theta2,real(f2)'),colorbar,caxis([-1 2])
    myplot2(['L = ' num2str(L)],'\phi (deg)','\theta (deg)')
        h5 = myaa(12);

end


%%


h6 = figure(3);clf,
scatter(180/pi*reshape(P1,[],1),180/pi*reshape(T1,[],1),'filled')
ylim([0 180])
xlim([0 360])
set(gca,'ydir','reverse')
myplot(['L = ' num2str(L)],'\phi (deg)','\theta (deg)')






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





%% vector filter


L = 8;
I = 2*L + 1;
J = L + 1;
phi = 2*pi*(0:(I-1))/I;
[muj, wj] = gqwt(J);
theta = acos(muj);
[T1 P1] = meshgrid(theta,phi);
tot1 = L^2 + 2*L;

[Bth Bphi Cth Cphi] = BC(L,T1,P1);

blm1 = zeros(tot1,1);
clm1 = zeros(tot1,1);
rng(8)
blm1 = rand(tot1,1) + 1i*rand(tot1,1);
clm1 = rand(tot1,1) + 1i*rand(tot1,1);

% blm1(7) = 1;
% clm1(1) = 0;

[Fth Fphi] = BCmult(blm1,clm1,Bth,Bphi,Cth,Cphi);
Fth = reshape(Fth,size(T1));
Fphi = reshape(Fphi,size(T1));

if 0
    figure(1),clf,imagesc(real(Fth)),colorbar;
    figure(2),clf,imagesc(imag(Fth)),colorbar;
    figure(3),clf,imagesc(real(Fphi)),colorbar;
    figure(4),clf,imagesc(imag(Fphi)),colorbar;
end


%%


[blm, clm] = vst(Fth,Fphi,L,muj,wj);

tmp = [blm1 blm];
figure(1),clf,plot([real(tmp)])
figure(2),clf,plot([imag(tmp)])
tmp = [clm1 clm];
figure(3),clf,plot([real(tmp)])
figure(4),clf,plot([imag(tmp)])

%%

[Fth2,Fphi2] = ivst(blm,clm,L,muj);



if 0
    figure(1),clf,imagesc(real(Fth-Fth2)),colorbar;
    figure(2),clf,imagesc(imag(Fth-Fth2)),colorbar;
    figure(3),clf,imagesc(real(Fphi-Fphi2)),colorbar;
    figure(4),clf,imagesc(imag(Fphi-Fphi2)),colorbar;
end


if 1
    figure(1),clf,imagesc(real(Fth)),colorbar;
    figure(2),clf,imagesc(imag(Fth)),colorbar;
    figure(3),clf,imagesc(real(Fphi)),colorbar;
    figure(4),clf,imagesc(imag(Fphi)),colorbar;
    figure(5),clf,imagesc(real(Fth2)),colorbar;
    figure(6),clf,imagesc(imag(Fth2)),colorbar;
    figure(7),clf,imagesc(real(Fphi2)),colorbar;
    figure(8),clf,imagesc(imag(Fphi2)),colorbar;
end

%%




%% fast vector filter basic


L = 20;
I = 2*L + 1;
J = L + 1;
phi = 2*pi*(0:(I-1))/I;
[muj, wj] = gqwt(J);
theta = acos(muj);
[T1 P1] = meshgrid(theta,phi);
tot1 = L^2 + 2*L;

[Bth Bphi Cth Cphi] = BC(L,T1,P1);

blm1 = zeros(tot1,1);
clm1 = zeros(tot1,1);
%rng(100)
%blm1 = -0.5 - 0.5*1i + rand(tot1,1) + 1i*rand(tot1,1);
%clm1 = -0.5 - 0.5*1i + rand(tot1,1) + 1i*rand(tot1,1);

blm1(lm2ind(K,3)) = 1;
clm1(lm2ind(K,3)) = 1;
%blm1(lm2ind(K-1,3)) = 1;
%clm1(lm2ind(K-1,3)) = 1;
%blm1(lm2ind(K+1,3)) = 1;
%clm1(lm2ind(K+1,3)) = 1;
%blm1(lm2ind(K-2,1)) = 1;
%clm1(lm2ind(K-2,1)) = 0;

[Fth Fphi] = BCmult(blm1,clm1,Bth,Bphi,Cth,Cphi);
Fth = reshape(Fth,size(T1));
Fphi = reshape(Fphi,size(T1));


K = 6;
P = 2*K + 1;
Q = K + 1;
phi = 2*pi*(0:(P-1))/P;
[muk, wk] = gqwt(Q);
theta = acos(muk);
[T2 P2] = meshgrid(theta,phi);
tot2 = K^2 + 2*K;

% truncatation method
blm2 = blm1(1:tot2);
clm2 = clm1(1:tot2);
[Fth2,Fphi2] = ivst(blm2,clm2,K,muk);

% fast scalar filter with correction
[Fth3, Fphi3] = vstfilterbasic(Fth,Fphi,L,muj,wj,K,muk,T1,T2,P2);


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
figure(13),clf,imagesc(real(Fth2-Fth3)),colorbar;
figure(14),clf,imagesc(imag(Fth2-Fth3)),colorbar;
figure(15),clf,imagesc(real(Fphi2-Fphi3)),colorbar;
figure(16),clf,imagesc(imag(Fphi2-Fphi3)),colorbar;






%% local interpolation


L = 7;
I = 2*L + 1;
J = L + 1;
phi1 = 2*pi*(0:(I-1))/I;
[muj, wj] = gqwt(J);
theta1 = acos(muj);
[T1 P1] = meshgrid(theta1,phi1);
tot1 = L^2 + 2*L;


blm1 = zeros(tot1,1);
clm1 = zeros(tot1,1);
%rng(100)
%blm1 = -0.5 - 0.5*1i + rand(tot1,1) + 1i*rand(tot1,1);
%clm1 = -0.5 - 0.5*1i + rand(tot1,1) + 1i*rand(tot1,1);

blm1(lm2ind(4,1)) = 0;
clm1(lm2ind(4,1)) = 1;
%blm1(lm2ind(L-1,3)) = 1;
%clm1(lm2ind(L-1,3)) = 1;
%blm1(lm2ind(K-2,1)) = 1;
%clm1(lm2ind(K-2,1)) = 0;

[Bth Bphi Cth Cphi] = BC(L,T1,P1);
[Fth Fphi] = BCmult(blm1,clm1,Bth,Bphi,Cth,Cphi);
Fth = reshape(Fth,size(T1));
Fphi = reshape(Fphi,size(T1));


K = 30;
P = 2*K + 1;
Q = K + 1;
phi2 = 2*pi*(0:(P-1))/P;
[muk, wk] = gqwt(Q);
theta2 = acos(muk);
[T2 P2] = meshgrid(theta2,phi2);
tot2 = K^2 + 2*K;

blm2 = zeros(tot2,1);
clm2 = zeros(tot2,1);
blm2(1:tot1) = blm1;
clm2(1:tot1) = clm1;

[Bth Bphi Cth Cphi] = BC(K,T2,P2);
[Fth2 Fphi2] = BCmult(blm2,clm2,Bth,Bphi,Cth,Cphi);
Fth2 = reshape(Fth2,size(T2));
Fphi2 = reshape(Fphi2,size(T2));

%%


L = 7;
I = 2*L + 1;
J = L + 1;
phi1 = 2*pi*(0:(I-1))/I;
[muj, wj] = gqwt(J);
theta1 = acos(muj);
[T1 P1] = meshgrid(theta1,phi1);
tot1 = L^2 + 2*L;

flm = zeros(tot1,1);
flm(lm2ind(7,3)) = 1;

ylm = sphericalY(L,T1,P1);
f1 = sum(repmat(flm,[1,numel(T1)]).*ylm,1);
f1 = reshape(f1,size(T1));



K = 30;
P = 2*K + 1;
Q = K + 1;
phi2 = 2*pi*(0:(P-1))/P;
[muk, wk] = gqwt(Q);
theta2 = acos(muk);
[T2 P2] = meshgrid(theta2,phi2);
tot2 = K^2 + 2*K;

flm2 = zeros(tot2,1);
flm2(1:tot1) = flm;

ylm = sphericalY(K,T2,P2);
f2 = sum(repmat(flm2,[1,numel(T2)]).*ylm,1);
f2 = reshape(f2,size(T2));


 



%%

%y1 = flipud(imag(Fth(1,:))');
%y2 = flipud(imag(Fth2(1,:))');
y1 = real(f1(1,:))'./sin(theta1);
y2 = real(f2(1,:))';

x1 = cos((theta1));
x2 = cos((theta2));
xtest = [-1; x1; 1];

figure(20),clf,hold all
plot(x1,y1,'o')
plot(x2,y2,'o')
hold off

N2 = length(x2);
p = 4;
yi = zeros(N2,1);
for m=1:length(x2),

xi = x2(m);

reg = find((sign(xtest-xi) > 0),1,'first')-1;

N1 = length(x1);
N2 = length(x2);


if N1 <= 4
    lowerlimit = 1;
    upperlimit = N1;
else
    if reg <= p
        lowerlimit = 1;
        upperlimit = 2*p;
    elseif reg > N1 - p
        lowerlimit = N1 - 2*p + 1;
        upperlimit = N1;
    else
        lowerlimit = reg - p;
        upperlimit = reg + p - 1;
    end
end


grab = lowerlimit:upperlimit;
u1 = zeros(length(grab),1);
for j=grab,
    tmp = 1;
    for k=grab,
        if k~=j
            tmp = tmp.*((xi-x1(k))/(x1(j)-x1(k)));
        end    
    end
    u1(j-grab(1) + 1) = tmp;
end

yi(m) = sum(y1(grab).*u1);


end

yi = yi.*sin(theta2);

figure(20),clf,hold all
plot(x1,y1,'o')
plot(x2,y2,'o')
plot(x2,yi,'o')
hold off

figure(21),clf,plot(yi-y2)



%%

figure(1),clf,imagesc(real(f1)),colorbar;
figure(2),clf,imagesc(imag(f1)),colorbar;
figure(3),clf,imagesc(real(f2)),colorbar;
figure(4),clf,imagesc(imag(f2)),colorbar;

figure(5),clf,imagesc(real(Fth2)),colorbar;
figure(6),clf,imagesc(imag(Fth2)),colorbar;
figure(7),clf,imagesc(real(Fphi2)),colorbar;
figure(8),clf,imagesc(imag(Fphi2)),colorbar;

%%



figure(1),clf,imagesc(real(Fth)),colorbar;
figure(2),clf,imagesc(imag(Fth)),colorbar;
figure(3),clf,imagesc(real(Fphi)),colorbar;
figure(4),clf,imagesc(imag(Fphi)),colorbar;

figure(5),clf,imagesc(real(Fth2)),colorbar;
figure(6),clf,imagesc(imag(Fth2)),colorbar;
figure(7),clf,imagesc(real(Fphi2)),colorbar;
figure(8),clf,imagesc(imag(Fphi2)),colorbar;

    
%%
    
    
    
    


tmp1 = Fth(1,:);

figure(1),clf,

%%







%%

p = 2;
theta2p = [-flipud(theta2(1:p)); theta2; pi + theta2(1:p)];

ts = pi;

ss = (sign(theta2p-ts) > 0);
ind = find(ss,1,'first')-1;

[ss theta2p];


Fth_pad = [fliplr(Fth(:,1:2)) Fth fliplr(Fth(:,(J-1):J))];
Fphi_pad = [fliplr(Fphi(:,1:2)) Fphi fliplr(Fphi(:,(J-1):J))];


figure(1),clf,imagesc(real(Fth)),colorbar;
figure(2),clf,imagesc(imag(Fth)),colorbar;
figure(3),clf,imagesc(real(Fphi)),colorbar;
figure(4),clf,imagesc(imag(Fphi)),colorbar;
figure(5),clf,imagesc(real(Fth2)),colorbar;
figure(6),clf,imagesc(imag(Fth2)),colorbar;
figure(7),clf,imagesc(real(Fphi2)),colorbar;
figure(8),clf,imagesc(imag(Fphi2)),colorbar;
figure(9),clf,imagesc(real(Fth_pad)),colorbar;
figure(10),clf,imagesc(imag(Fth_pad)),colorbar;
figure(11),clf,imagesc(real(Fphi_pad)),colorbar;
figure(12),clf,imagesc(imag(Fphi_pad)),colorbar;







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
f3 = fssbasic(f,Level(l+1).L,Level(l).L);
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
f3 = fssbasic(f,Level(n+1).L,Level(n).L);
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



