



%% smatrix_thin_circular_cylinder

k = 2*pi;
a = 0.01;
L = 0.3;
er = 2;

th_i = linspace(0,pi,20);
phi_i = linspace(0,2*pi,20);
[Ti Pi] = ndgrid(th_i,phi_i);

th_s = linspace(0,pi,20);
phi_s = linspace(0,2*pi,20);
[Ts Ps] = ndgrid(th_s,phi_s);

[Svv Svh Shv Shh] = smatrix_thin_circular_cylinder(k,a,L,er,Ts,Ps,Ti,Pi);

imagesc(Svv(:,:,1,1)),colorbar
imagesc(Svh(:,:,1,1)),colorbar
imagesc(Shv(:,:,1,1)),colorbar
imagesc(Shh(:,:,1,1)),colorbar







