
%% momGmatrix2D

lam = 1;
k = 2*pi/lam;

% create a dense square grid of points
N = 40;
x = linspace(-2,2,N);
y = linspace(-2,2,N);
[X Y] = meshgrid(x,y);
dx = x(2)-x(1);


R = sqrt(X.^2 + Y.^2);
o = k^2*(1+1i*0.2);
C = zeros(N);

rng(1);
I = randi(N^2,[N^2/2 1]);
C(I(:)) = o;

figure(1),clf
imagesc(x,y,real(C)),colorbar
set(gca,'ydir','normal')

% incident field
sx = 4;
sy = 4;

ind = find(C > 0);
O = C(ind);
NN = length(ind);
X2 = X(ind);
Y2 = Y(ind);

% return full matrix, set primed and unprimed coordinates to be the same
Xp = X2;
Yp = Y2;

% build the matrix
G = momGmatrix2D(dx,k,X2,Y2,Xp,Yp);

% incident field 
phi_inc = exp(1i*k*sqrt((X2-sx).^2 + (Y2-sy).^2));

% mom solution
phi = (eye(NN) - G*diag(O))\phi_inc;

% remap
Phi = zeros(size(C));
Phi(ind) = phi;

% plot
figure(2),clf
imagesc(x,y,real(Phi)),colorbar
set(gca,'ydir','normal')

% scattered field 
sca = sum(phi.*O.*phi_inc)



%%

[Q D] = eig(G,'vector');

u = Q.'*phi_inc;
v = Q\phi_inc;
w = o./(1-o*D);

dji = sum(u.*w.*v)




%%

[Q D] = eig(G);

imagesc(abs(Q*D/Q))

N = length(G(:,1));

C = ones(N/4,1);
t = 1./[(2 + 2*1i)*C; (3.5 + 1i)*C; (2.5 + 0.2*1i)*C; (4 + 0.01*1i)*C];
T = diag(t);

B = inv(Q)*T*Q;

imagesc(abs(B)),colorbar

plot(diag(B))
mean(diag(B))

C = inv(B + D);

imagesc(abs(C)),colorbar

plot([real(diag(C)) imag(diag(C))])


% 
% 
% [U S V] = svd(B);
% plot(diag(S))
% 
% grab = 100;
% F = U(:,1:grab)*S(1:grab,1:grab)*V(:,1:grab)';
% imagesc(abs(F))





%%

[U S V] = svd(G);

imagesc(abs(G))
imagesc(abs(U*S*V')),colorbar

A = 3*randn(size(G)) + 1i*randn(size(G));

B = V*A*U';
imagesc(abs(B)),colorbar





