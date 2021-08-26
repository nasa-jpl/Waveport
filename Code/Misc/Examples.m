
%% fftfreq


% even number of points
fs = 1;
N = 10;
f = fftfreq(fs,N)

figure(1),clf
stem(f)
myplot('fftfreq','DFT Index','f (rad)')
ylim([-fs/2 fs/2])

% odd number of points
fs = 10;
N = 11;
f = fftfreq(fs,N)

figure(1),clf
stem(f)
myplot('fftfreq','DFT Index','f (rad)')
ylim([-fs/2 fs/2])


%% phase wrapping

x = linspace(0,20,1000);
phi = 2*pi*x;

phase1 = angle(exp(1i*f));
phase2 = wrap(f);

mean(abs(phase1-phase2))


    

%% cubicroots

%%% roots for multiple polynomials
Npolynomials = 3;
rng(2);
coef = rand(Npolynomials,4);  % generate coefficients

% compute vectorized roots using analytic solution
[x1 x2 x3] = cubicroots(coef(:,1),coef(:,2),coef(:,3),coef(:,4));
wp_roots = [x1 x2 x3];

%%% compute roots with matlab's roots function
matlab_roots = zeros(Npolynomials,3);
for n=1:Npolynomials;
    matlab_roots(n,:) = roots(coef(n,:));    
end

%%% they are in different arraigments, so just plot all roots on top of each other
figure(1),clf
hold all
plot(real(wp_roots(:)),imag(wp_roots(:)),'o')
plot(real(matlab_roots(:)),imag(matlab_roots(:)),'.')
hold off
myplot('Cubic roots','Real','Imag')



%% quarticroots

%%% roots for multiple polynomials
Npolynomials = 3;
rng(3);
coef = rand(Npolynomials,5);  % generate coefficients

% compute vectorized roots using analytic solution
[x1 x2 x3 x4] = quarticroots(coef(:,1),coef(:,2),coef(:,3),coef(:,4),coef(:,5));
wp_roots = [x1 x2 x3 x4];

%%% compute roots with matlab's roots function
matlab_roots = zeros(Npolynomials,4);
for n=1:Npolynomials;
    matlab_roots(n,:) = roots(coef(n,:));    
end

%%% they are in different arraigments, so just plot all roots on top of each other
figure(1),clf
hold all
plot(real(wp_roots(:)),imag(wp_roots(:)),'o')
plot(real(matlab_roots(:)),imag(matlab_roots(:)),'.')
hold off
myplot('Quartic roots','Real','Imag')



%% discoball spherical sampling

Nlat = 20;
[th phi] = discoball(Nlat);
[x y z] = sph2cart(ones(size(th)),th,phi);

h=figure(1);,clf
scatter3(x,y,z,40,z,'filled')
axis equal
myplot('Disco ball','x','y','z')
view([1 -1 0.4])



%% random spherical points

rng(10)
N = 1000;
[th phi] = randsphere(N);
[x y z] = sph2cart(ones(size(th)),th,phi);

h=figure(1);,clf
scatter3(x,y,z,40,z,'filled')
axis equal
myplot('Uniformily Random Points on a Sphere','x','y','z')
view([1 -1 0.4])



