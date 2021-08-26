




%%  Test scalar axial translation, alphaz

% set the background wave number equal to 1
k = 1;

% create a local cartesian grid for frame j
Np = 10;
ext = 2;
x = linspace(-ext,ext,Np);
y = x;
z = x;
[Xj Yj Zj] = ndgrid(x,y,z);

% magnitude of z translation from frame i to j
rji = 100;

% coordinates of the same points in frame i
Xi = Xj;
Yi = Yj;
Zi = Zj + rji;

% convert the points to spherical coordinates in each frame
[Ri Ti Pi] = cart2sph(Xi,Yi,Zi);
[Rj Tj Pj] = cart2sph(Xj,Yj,Zj);

% maximum degree harmonic in the translation matrix, L and L'
L = 15;
Lp = 15;

% evalute the scalar wave functions at those points relative to each frame
% outgoing waves from frame i
% incoming waves in frame j (regular)
psi_i = psilm(Lp,k,Ri,Ti,Pi);
psi_j = psilm(L,k,Rj,Tj,Pj,'rg');

% total number of harmonics
tot = Lp^2+2*Lp+1;

% create expansion coefficients for the outgoing waves in frame i
rng(1);
alm = randn(tot,1) + 1i*randn(tot,1);

% compute the axial translation matrix
alp = alphaz(L,Lp,k,rji);

if 0
    % plot the axial translation matrix
    figure(1),clf
    imagesc(db(alp)),colorbar
    myplot('Scalar Axial Translation Matirx','(l'',m'') Linear Index','(l,m) Linear Index')
    mycbar('(dB)',[-200 -20])
end

% evaluate the incoming expansion coefficients in frame j
blm = alp*alm;

% evaluate the expansions of the field in each frame
phi_i = psi_i*alm;
phi_j = psi_j*blm;

% compute the average error
mean(abs(phi_i-phi_j))

% plot the comparison
% change the number of L and Lp to see the effect of the number of
% harmonics on the accuracy of the translation
figure(2),clf,hold all
plot([real(phi_i) imag(phi_i)])
plot([real(phi_j) imag(phi_j)],'--')
hold off
myplot('Scalar field','Grid point','Field')
legend('real(\phi_i)','imag(\phi_i)','real(\phi_j)','imag(\phi_j)')

% plot the residual error
figure(3),clf
plot([real(phi_i-phi_j) imag(phi_i-phi_j)])
myplot('Residual Error','Grid point','\phi_i-\phi_j')
legend('Real','Imag')



%% Test sparse scalar axial translation, alphazSparse

% set wavenumber and translation distance
k = 1;
rji = 100;

% maximum degree harmonic in the translation in each frame
L = 21;
Lp = 9;

% computed on full matrix
alp = alphaz(L,Lp,k,rji);

% sparse matrix entries for the same
[row col alpsp] = alphazSparse(L,Lp,k,rji);
[m n] = size(alp); 
% create sparse matrix
alp2 = sparse(row,col,alpsp,m,n);

% compare
mean(abs(alp(:) - alp2(:)))



%% Test scalar translation matrix, scalarTranslation

% set wavenumber
k = 1;

% create grid of points for frame j
Np = 10;
ext = 2;
x = linspace(-ext,ext,Np);
y = x;
z = x;
[Xj Yj Zj] = ndgrid(x,y,z);

% define the translation vector that points from frame i to frame j
rji = [-100 100 -100]

% coordinates of the same points in frame i
Xi = Xj + rji(1);
Yi = Yj + rji(2);
Zi = Zj + rji(3);

% convert to spherical coordinates
[Ri Ti Pi] = cart2sph(Xi,Yi,Zi);
[Rj Tj Pj] = cart2sph(Xj,Yj,Zj);

% maximum degree harmonics for the translation, L and L'
L = 10;
Lp = 10;

% evalute the scalar wave functions at those points relative to each frame
% outgoing waves from frame i
% incoming waves in frame j (regular)
psi_i = psilm(Lp,k,Ri,Ti,Pi);
psi_j = psilm(L,k,Rj,Tj,Pj,'rg');

% total number of harmonics in frame i
tot = Lp^2+2*Lp+1;

% create expansion coefficients for the outgoing waves in frame i
rng(1);
alm = randn(tot,1) + 1i*randn(tot,1);

% apply translation matrix to the coefficients
[blm] = scalarTranslation(alm,L,Lp,k,rji,[]);

% compute fields
phi_i = psi_i*alm;
phi_j = psi_j*blm;

% compute the average error
mean(abs(phi_i-phi_j))

% plot the comparison
% change the number of L and Lp to see the effect of the number of
% harmonics on the accuracy of the translation
figure(2),clf,hold all
plot([real(phi_i) imag(phi_i)])
plot([real(phi_j) imag(phi_j)],'--')
hold off
myplot('Scalar field','Grid point','Field')
legend('real(\phi_i)','imag(\phi_i)','real(\phi_j)','imag(\phi_j)')

% plot the residual error
figure(3),clf
plot([real(phi_i-phi_j) imag(phi_i-phi_j)])
myplot('Residual Error','Grid point','\phi_i-\phi_j')
legend('Real','Imag')



%% Test full scalar translation matrix, alpha

% unit background wavenumber
k = 1;

% arbitrary direction
rji = 20*[-1 0.5 -0.75];

% a few different tests
for test = 1:3,  % loop over different matrix sizes
for flavor = 1:2; % loop over 'rg' and not 'rg'
    switch test
        case 1
            L = 5;
            Lp = 10;
        case 2
            L = 10;
            Lp = 10;
        case 3
            L = 10;
            Lp = 5;
    end
    switch flavor
        case 1
            str = [];
        case 2
            str = 'rg';
    end

    % create blank coefficient matrix
    tot = Lp^2 + 2*Lp + 1;
    alm = ones(tot,1);

    % compute alpha using the diagonalized form
    [~, alp1] = scalarTranslation(alm,L,Lp,k,rji,str);

    % compute alpha from the full matrix recursion
    alp2 = alpha(L,Lp,k,rji,str);

    % check the comparision
    mean(abs(alp1(:) - alp2(:)))
end
end

figure(1),clf
imagesc(real(alp1)),colorbar
myplot('Scalar Matrix: Diagonalized','Linear (l'',m'')','Linear (l,m)')

figure(2),clf
imagesc(real(alp2)),colorbar
myplot('Scalar Matrix: Full','Linear (l'',m'')','Linear (l,m)')

figure(3),clf
imagesc(real(alp1-alp2)),colorbar
myplot('Scalar Matrix: Diagonalized - Full','Linear (l'',m'')','Linear (l,m)')



%% Test vector axial translation, AzBz

% set the background wave number equal to 1
k = 1;

% create a local cartesian grid for frame j
Np = 10;
ext = 2;
x = linspace(-ext,ext,Np);
y = x;
z = x;
[Xj Yj Zj] = ndgrid(x,y,z);

% magnitude of z translation from frame i to j
rji = 100;

% coordinates of the same points in frame i
Xi = Xj;
Yi = Yj;
Zi = Zj + rji;

% convert the points to spherical coordinates in each frame
[Ri Ti Pi] = cart2sph(Xi,Yi,Zi);
[Rj Tj Pj] = cart2sph(Xj,Yj,Zj);

% maximum degree harmonic in the translation matrix, L and L'
L = 20;
Lp = 25;

% set if use normalized vector wave functions
normstr = [];

% evaluate vector spherical wave functions in frame i and j
[Mth_i Mphi_i Nr_i Nth_i Nphi_i] = MN(Lp,k,Ri,Ti,Pi,[],[],normstr);
[Mth_j Mphi_j Nr_j Nth_j Nphi_j] = MN(L,k,Rj,Tj,Pj,'rg',[],normstr);

% total number of harmonics in frame i
tot = Lp^2+2*Lp;

% create expansion coefficients for the outgoing waves in frame i
rng(1);
alm = randn(tot,1) + 1i*randn(tot,1);
blm = randn(tot,1) + 1i*randn(tot,1);

% evaluate the axial translation matrix
[A B] = AzBz(L,Lp,k,rji,[],normstr);

% compute the incoming expansion coefficients in frame j
clm = A*alm + B*blm;
dlm = B*alm + A*blm;

% evaluate the vector field in each frame
[Er_i Eth_i Ephi_i] = MNmult(Mth_i,Mphi_i,Nr_i,Nth_i,Nphi_i,alm,blm);
[Er_j Eth_j Ephi_j] = MNmult(Mth_j,Mphi_j,Nr_j,Nth_j,Nphi_j,clm,dlm);

% convert the vector components to Cartesian for comparison
[Exi Eyi Ezi] = sph2cart(Ri(:),Ti(:),Pi(:),Er_i,Eth_i,Ephi_i);
[Exj Eyj Ezj] = sph2cart(Rj(:),Tj(:),Pj(:),Er_j,Eth_j,Ephi_j);

if 0
    % plot the vector field evaluated in frame i
    figure(1),clf,
    quiver3(Xi(:),Yi(:),Zi(:),real(Exi),real(Eyi),real(Ezi));
    myplot('Field in Frame i','x','y','z')

    % plot the vector field evaluated in frame j
    figure(2),clf,
    quiver3(Xj(:),Yj(:),Zj(:),real(Exj),real(Eyj),real(Ezj));
    myplot('Field in Frame j','x','y','z')
end

% compare
mean(abs(Exi-Exj))
mean(abs(Eyi-Eyj))
mean(abs(Ezi-Ezj))


%% Test sparse vector axial translation, AzBzSparse

% set wavenumber and translation distance
k = 1;
rji = 100;

% maximum degree harmonic in the translation in each frame
L = 21;
Lp = 13;

% set normalization and regular strings
normstr = [];
rgstr = [];

% computed on full matrix
[A B] = AzBz(L,Lp,k,rji,rgstr,normstr);

% % sparse matrix entries for the same
[row col A2 B2] = AzBzSparse(L,Lp,k,rji,rgstr,normstr);
[m n] = size(A);

% convert to sparse matrices
A2 = sparse(row,col,A2,m,n);
B2 = sparse(row,col,B2,m,n);

% compuare
mean(abs(A(:) - A2(:)))
mean(abs(B(:) - B2(:)))




%% Test vector translation

% set the background wave number equal to 1
k = 1;

% create a local cartesian grid for frame j
Np = 10;
ext = 2;
x = linspace(-ext,ext,Np);
y = x;
z = x;
[Xj Yj Zj] = ndgrid(x,y,z);

% define the translation vector that points from frame i to frame j
rji = [-100 100 -100]

% coordinates of the same points in frame i
Xi = Xj + rji(1);
Yi = Yj + rji(2);
Zi = Zj + rji(3);

% convert the points to spherical coordinates in each frame
[Ri Ti Pi] = cart2sph(Xi,Yi,Zi);
[Rj Tj Pj] = cart2sph(Xj,Yj,Zj);

% maximum degree harmonic in the translation matrix, L and L'
L = 20;
Lp = 25;

% set if use normalized vector wave functions
normstr = [];

% evaluate vector spherical wave functions in frame i and j
[Mth_i Mphi_i Nr_i Nth_i Nphi_i] = MN(Lp,k,Ri,Ti,Pi,[],[],normstr);
[Mth_j Mphi_j Nr_j Nth_j Nphi_j] = MN(L,k,Rj,Tj,Pj,'rg',[],normstr);

% total number of harmonics in frame i
tot = Lp^2+2*Lp;

% create expansion coefficients for the outgoing waves in frame i
rng(1);
alm = randn(tot,1) + 1i*randn(tot,1);
blm = randn(tot,1) + 1i*randn(tot,1);

% apply translation matrix to the coefficients
[clm dlm] = vectorTranslation(alm,blm,L,Lp,k,rji,[],normstr);

% evaluate the vector field in each frame
[Er_i Eth_i Ephi_i] = MNmult(Mth_i,Mphi_i,Nr_i,Nth_i,Nphi_i,alm,blm);
[Er_j Eth_j Ephi_j] = MNmult(Mth_j,Mphi_j,Nr_j,Nth_j,Nphi_j,clm,dlm);

% convert the vector components to Cartesian for comparison
[Exi Eyi Ezi] = sph2cart(Ri(:),Ti(:),Pi(:),Er_i,Eth_i,Ephi_i);
[Exj Eyj Ezj] = sph2cart(Rj(:),Tj(:),Pj(:),Er_j,Eth_j,Ephi_j);

if 1
    % plot the vector field evaluated in frame i
    figure(1),clf,
    quiver3(Xi(:),Yi(:),Zi(:),real(Exi),real(Eyi),real(Ezi));
    myplot('Field in Frame i','x','y','z')

    % plot the vector field evaluated in frame j
    figure(2),clf,
    quiver3(Xj(:),Yj(:),Zj(:),real(Exj),real(Eyj),real(Ezj));
    myplot('Field in Frame j','x','y','z')
end

% compare
mean(abs(Exi-Exj))
mean(abs(Eyi-Eyj))
mean(abs(Ezi-Ezj))


%% Test the full vector translation matrix, AB

% unit background wavenumber
k = 1;

% arbitrary direction
rji = 20*[1 -0.5 -0.75];

% a few different tests
for test = 1:3,     % loop over different matrix sizes
for flavor = 1:2;   % loop over 'rg' and not rg
for spin = 1:2;    % loop over 'norm' and not norm
    switch test
        case 1
            L = 8;
            Lp = 10;
        case 2
            L = 10;
            Lp = 10;
        case 3
            L = 10;
            Lp = 8;
    end
    switch flavor
        case 1
            rgstr = [];
        case 2
            rgstr = 'rg';
    end
    switch spin
        case 1
            normstr = [];
        case 2
            normstr = 'norm';
    end
        
    % create blank coefficient matrix
    tot = Lp^2 + 2*Lp;
    alm = ones(tot,1);
    blm = ones(tot,1);

    % compute A and B using the diagonalized form
    [~,~,A1,B1] = vectorTranslation(alm,blm,L,Lp,k,rji,rgstr,normstr);

    % compute A and B from the full matrix recursion
    [A2 B2] = AB(L,Lp,k,rji,rgstr,normstr);

    % check the comparision
    [mean(abs(A1(:) - A2(:))) mean(abs(B1(:) - B2(:)))]

end
end
end

figure(1),clf
imagesc(real(A1)),colorbar
myplot('Vector Matrix: Diagonalized','Linear (l'',m'')','Linear (l,m)')

figure(2),clf
imagesc(real(A2)),colorbar
myplot('Vector Matrix: Full','Linear (l'',m'')','Linear (l,m)')

figure(3),clf
imagesc(real(A1-A2)),colorbar
myplot('Vector Matrix: Diagonalized - Full','Linear (l'',m'')','Linear (l,m)')





