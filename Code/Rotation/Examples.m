

%% Euler angle conversions 

% Euler angles
alpha = pi/6;
beta = pi/5;
gamma = pi/4;

% ZXZ rotation matrix from Euler angles
R = euler2rot(alpha,beta,gamma);

% Euler angles from rotation matrix
[alpha2 beta2 gamma2] = rot2euler(R);


%% ZXZ Rotation

% Euler angles
alpha = pi/6;
beta = pi/5;
gamma = pi/4;

% successive rotations
R0 = euler2rot(0,0,0);
R1 = euler2rot(alpha,0,0);
R2 = euler2rot(alpha,beta,0);
R3 = euler2rot(alpha,beta,gamma);

% plot the origal frame (black) and rotated frame (rbg)
vw = [1 1 0.8];
sz = 1;
origin = [0 0 0];
figure(1),clf
plotrot(R0,origin,sz,'k',':')
plotrot(R1,origin,sz,[],'--')
view(vw)
axis equal
axis([-1 1 -1 1 -1 1])
myplot('Z Rotation','x','y','z')

figure(1),clf
plotrot(R1,origin,sz,'k',':')
plotrot(R2,origin,sz,[],'--')
view(vw)
axis equal
axis([-1 1 -1 1 -1 1])
myplot('X Rotation','x','y','z')

figure(1),clf
plotrot(R2,origin,sz,'k',':')
plotrot(R3,origin,sz,[],'--')
view(vw)
axis equal
axis([-1 1 -1 1 -1 1])
myplot('Z Rotation','x','y','z')


%% Dlmp

% create rotation
L = 6;
alpha = pi/2;
beta = pi/2;
gamma = pi/3;

% compute Dlmp (inverse rotation) up to harmonics L
D = Dlmp(L,alpha,beta,gamma);

figure(1),clf
imagesc(abs(D))
myplot('D_{lmp}(\alpha,\beta,\gamma)','(l,p) linear index','(l,m) linear index')




%% Test Dlmp by rotating spherical harmoncis

% create spherical grid of points
N = 30;
theta = linspace(0,pi,N);
phi = linspace(0,2*pi,2*N);
[T P] = ndgrid(theta,phi);

% convert to cartesian
[X Y Z] = sph2cart(ones(size(T)),T,P);

% choose maximum degree
L = 6;
tot = L^2 + 2*L;

% compute spherical harmonics at all points
ylm = sphericalY(L,T,P);

% create arbitrary expansion coefficients
rng(1);
alm = randn(tot,1) + 1i*rand(tot,1);

% compute the spherical function from the coefficients
f = reshape(ylm*alm,size(T));

% plot the function
figure(1),clf
surf(X,Y,Z,real(f))
shading flat
axis(1.5*[-1 1 -1 1 -1 1])
hold on
plotrot(diag([1 1 1]),[0 0 0],1.5)
hold off
view([1 1 1])
myplot('Spherical Function','x','y','z');


% create rotation
alpha = pi/5;
beta = pi/3;
gamma = 0;

% compute Dlmp (inverse rotation) up to harmonics L
D = Dlmp(L,alpha,beta,gamma);



% choose forward or inverse and compute associated rotation matrix
if 1 % forward
    D = D'; % need to conjugate the output of Dlmp
    R = euler2rot(alpha,beta,gamma);
else % inverse
	R = euler2rot(alpha,beta,gamma)';
end

% rotation the expansion coefficents
almp = D*alm;

% compute spherical function from rotated coefficients
f2 = reshape(ylm*almp,size(T));

% plot the rotated function
figure(2),clf
surf(X,Y,Z,real(f2))
shading flat
axis(1.5*[-1 1 -1 1 -1 1])
hold on
plotrot(R,[0 0 0],1.5)
hold off
view([1 1 1])
myplot('Spherical Function - Rotated','x','y','z');



%% Test DlmpSparse

L = 8;
alpha = pi/7;
beta = pi/8;
gamma = pi/3;

% compute Dlmp
D = Dlmp(L,alpha,beta,gamma);

% compute sparse Dlmp
[row col DS] = DlmpSparse(L,alpha,beta,gamma);

% low sparse matrix into a full matrix
D2 = zeros(size(D));
D2(sub2ind(size(D),row,col)) = DS(:);

% compare the two
mean(abs(D(:) - D2(:)))


%% Compare to DlmpSparseFast 

% pick a large L
L = 100;
alpha = pi/3;
beta = pi/4;
gamma = pi/5;

% compute DlmpSparse and time
tic
[row col Dlmps] = DlmpSparse(L,alpha,beta,gamma);
toc

% compute DlmpSparseFast and time
tic
[row2 col2 Dlmps2] = DlmpSparseFast(L,alpha,beta,gamma);
toc

% check that the results are the same
c1 = mean(abs(row2-row));
c2 = mean(abs(col2-col));
c3 = mean(abs(Dlmps2-Dlmps));
[c1 c2 c3]


%% dlmp, compare all the routines for \alpha = 0, \gamma = 0

% maximum degree harmonic and angles
L = 10;
beta = pi/3;

% compute the matrix directly
dlmp1 = dlmpBetaDirect(L,beta);

% compute it using the recursive algorithm
dlmp2 = dlmpBeta(L,beta);

% compute it using the faster version of the recursive algorithm
dlmp3 = dlmpBetaFast(L,beta);

% compute it using the sparse algorithm
[row col dlmpS] = dlmpBetaSparse(L,beta);

% load sparse matrix into full matrix for easy comparison
dlmp4 = zeros(size(dlmp1));
dlmp4(sub2ind(size(dlmp1),row,col)) = dlmpS(:);

% compute the full Dlmp with \alpha = 0, \gamma = 0
D = Dlmp(L,0,beta,0);

% compare them to the full Dlmp
mean(abs(dlmp1(:) - D(:)))
mean(abs(dlmp2(:) - D(:)))
mean(abs(dlmp3(:) - D(:)))
mean(abs(dlmp4(:) - D(:)))

% compare to each other
mean(abs(dlmp2(:) - dlmp3(:)))
mean(abs(dlmp2(:) - dlmp4(:)))


%% dlmp, Test 2, compare with full rotation matrix 

% maximum degree harmonic and angles
L = 10;
alpha = pi/3;
beta = pi/2;
gamma = pi/3;
tab = lmtable(L);
m = tab(:,2);

% compute full Dlmp
D = Dlmp(L,alpha,beta,gamma);

% compute little-d only using the fast full matrix version
dlmp5 = dlmpBetaFast(L,beta);

% multiply by the diagonal matricies for \alpha and \gamma rotations
D2 = diag(exp(1i*alpha*m))*dlmp5*diag(exp(1i*gamma*m));

% compare
mean(abs(D2(:)-D(:)))



%% Test Dlmp properties

% create rotation
L = 6;
alpha = pi/5;
beta = pi/3;
gamma = pi/3;

%%% test transpose relation 
D1 = Dlmp(L,alpha,beta,gamma);
D2 = Dlmp(L,-gamma+pi,beta,-alpha+pi)';

% compare
mean(abs(D1(:)-D2(:)))


%%% test the conjugation relation 
L = 6;
alpha = pi/5;
beta = pi/3;
gamma = pi/3;

D1 = Dlmp(L,alpha,beta,gamma);
D2 = conj(D1);

% create indices
tot = L^2 + 2*L;
[l m] = ind2lm(1:tot);
ind_minus_m = lm2ind(l,-m);
[M P] = ndgrid(m,m);

% apply constants to conjugated matrix
D2 = (-1).^(M+P).*D2;

% evaluate unconjugate matrix at negative indices
D1 = D1(:,ind_minus_m);
D1 = D1(ind_minus_m,:);

% compare
mean(abs(D1(:)-D2(:)))


%%% test Legendre relation
L = 6;
beta = pi/3;
d1 = dlmpBeta(L,beta);

% create indices for (l00)
l = 1:L;
ind = lm2ind(l,zeros(size(l)));
ind2 = sub2ind(size(d1),ind,ind);
% grab (l00) matrix elements
d1 = d1(ind2);

% compute Legendre polynomials
leg = legendrePl(L,cos(beta));
% exclude monopole
leg = leg(2:end);

% compare
mean(abs(d1(:)-leg(:)))


%%% test Dirichlet kernel
L = 6;
beta = pi/3;
d1 = dlmpBeta(L,beta);
d1 = diag(d1);

% compute sum of the diagonals of the blocks
diagsum = zeros(L,1);
for l=1:L,
    ind = lm2ind(l*ones(2*l+1,1),(-l:l)');
    diagsum(l) = sum(d1(ind));
end

% compute Dirichlet kernel
dirich = sin((2*(1:L)'+1)/2*beta)/sin(beta/2);

mean(abs(diagsum(:)-dirich(:)))



