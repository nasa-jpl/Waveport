
% Examples for Legendre polynomial routines


%% legendrePl

N = 100;
L = 5;
x = linspace(-1,1,N);
Pl = legendrePl(L,x);

% check against Matlab's legendre for m = 0
leg = zeros(L+1,N);
for l = 0:L,
    tmp = legendre(l,x);
    % grab m = 0
    leg(l+1,:) = tmp(1,:);
end

figure(1),clf,hold all
plot(x,Pl)
plot(x,leg','--')
hold off
myplot('Legendre P_l(x)','x',' ')



%% legendrePlp

N = 200;
L = 5;
x = linspace(-1,1,N);
Plp = legendrePlp(L,x);

% check against numeric gradient
dx = x(2)-x(1);
Pl = legendrePl(L,x);
Plp_numeric = gradient(Pl,dx);

figure(1),clf,hold all
plot(x,Plp')
plot(x,Plp_numeric','--')
hold off
myplot('Legendre P_l(x)'' ','x',' ')


%% \tilde Plm

N = 200;
L = 3;
x = linspace(-1,1,N);
plm = Plm(L,x);

figure(1),clf,hold all
plot(x,plm')
hold off
myplot('P^l_m(x)','x',' ')
tab = lmtable(L,'mono');
leg = legend(num2str(tab),'location','eastoutside')
myleg(leg,'       l   m')



%% d/dx \tilde Plm(x)

N = 500;
L = 3;
x = linspace(-1,1,N);
plmp = Plmp(L,x);

% check against numerical gradient
dx = x(2)-x(1);
plm = Plm(L,x);
plmp_numeric = gradient(plm,dx);

figure(1),clf,hold all
plot(x,plmp')
plot(x,plmp_numeric','--')
hold off
myplot('P_l^m(x)'' ','x',' ')
tab = lmtable(L,'mono');
leg = legend(num2str(tab),'location','eastoutside');
myleg(leg,'       l   m');


%% m \tilde Plm(cos(theta))/sin(theta)

N = 100;
L = 3;
theta = linspace(0+pi/200,pi-pi/200,N); % avoid the poles
mPsin = mPlmsin(L,theta);

% check against direction computation with Plm
Ltot = L^2 + 2*L;
mm = zeros(Ltot,1);
cnt = 1;
for l=1:L,
for m=-l:l,
    mm(cnt) = m;
    cnt = cnt+1;
end
end
[MM T] = ndgrid(mm,theta);
Plmtmp = Plm(L,cos(theta));
mPsin2 = MM./sin(T).*Plmtmp(2:end,:);

% plot
figure(1),clf,hold all
plot(theta,mPsin')
plot(theta,mPsin2','--')
hold off
myplot('m P_l^m(cos\theta)/sin\theta ','\theta (rad)',' ')
tab = lmtable(L);
leg = legend(num2str(tab),'location','eastoutside');
myleg(leg,'       l   m');



%% d/dtheta tilde Plm(cos(theta))

N = 500;
L = 3;
theta = linspace(0,pi,N);
dtheta_Plm = Plmp2(L,theta);

% check against numerical gradient
dt = theta(2)-theta(1);
plm = Plm(L,cos(theta));

% trim off monopole
plm = plm(2:end,:);
plmp_numeric = gradient(plm,dt);

figure(1),clf,hold all
plot(x,dtheta_Plm')
plot(x,plmp_numeric','--')
hold off
myplot('d/d\theta P_l^m(cos\theta)','\theta (rad)',' ')
tab = lmtable(L);
leg = legend(num2str(tab),'location','eastoutside');
myleg(leg,'       l   m');

