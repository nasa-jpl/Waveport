
%% spherical harmonic integration

L = 6;

I = L+1;
J = ceil((L+1)/2);
[xj wj] = legpts(J);
theta_j = acos(xj);
phi_i = (0:(I-1))*2*pi/I;
[P T] = meshgrid(phi_i,theta_j);
[~, W] = meshgrid(phi_i,wj);

ylm = sphericalY(L,T,P,'mono');
tot = L^2 + 2*L + 1;

quadint = 2*pi/(L+1)*ylm.*repmat(W(:),1,tot);
qi = sum(quadint,1);


J = 500;
I = 2*J;
theta = linspace(0,pi,J);
dtheta = theta(2)-theta(1);
dphi = 2*pi/I;
phi = dphi*(0:(I-1));
[T P] = ndgrid(theta,phi);

ylm = sphericalY(L,T,P,'mono');

q = zeros(1,tot);
for n=1:tot,
    tmp = reshape(ylm(:,n),size(T));
    tmp = trapz(theta,tmp.*sin(T),1);
    tmp = trapz(phi,tmp);
    q(n) = tmp;
end

plot(abs([qi - q]),'o')



%% integration of products of spherical harmonics 

L = 6;

I = 2*L+1;
J = L+1;
[xj wj] = legpts(J);
theta_j = acos(xj);
phi_i = (0:(I-1))*2*pi/I;
[P T] = meshgrid(phi_i,theta_j);
[~, W] = meshgrid(phi_i,wj);

ylm = sphericalY(L,T,P,'mono');
tot = L^2 + 2*L + 1;


quadint = 2*pi/I*ylm.^2.*repmat(W(:),1,tot);
qi = sum(quadint,1);


J = 600;
I = 2*J;
theta = linspace(0,pi,J);
dtheta = theta(2)-theta(1);
dphi = 2*pi/I;
phi = dphi*(0:(I-1));
[T P] = ndgrid(theta,phi);
ylm = sphericalY(L,T,P,'mono');

q = zeros(1,tot);
for n=1:tot,
    tmp = reshape(ylm(:,n).^2,size(T));
    tmp = trapz(theta,tmp.*sin(T),1);
    tmp = trapz(phi,tmp);
    q(n) = tmp;
end


plot(abs([qi' q']),'o')

plot(abs([qi - q]),'o')


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

%%


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


