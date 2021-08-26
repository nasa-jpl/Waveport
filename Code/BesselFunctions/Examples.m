
% Spherical Bessel function examples

%% jinc

x = linspace(0,40,1000)';
j1 = jinc(x);

figure(1),clf,hold all
plot(x,j1)
hold off
myplot('jinc(x)','x',' ')


%% sbesselj

x = linspace(0,20,1000)';
n = 0:5;
[X N] = ndgrid(x,n);
j1 = sbesselj(N,X);

figure(1),clf,hold all
plot(x,j1)
hold off
myplot('sbesselj(x)','x',' ')


%% sbesselj2

x = linspace(0,20,1000)';
n = 0:5;
[X N] = ndgrid(x,n);
j1 = sbesselj2(N,X);

figure(1),clf,hold all
plot(x,j1)
hold off
myplot('sbesselj2(x)','x',' ')
ylim([-0.1 0.5])


%% sbesseljp

x = linspace(0,20,1000)';
n = 0:5;
[X N] = ndgrid(x,n);
j1 = sbesseljp(N,X);

figure(1),clf,hold all
plot(x,j1)
hold off
myplot('sbesseljp(x)','x',' ')
ylim([-1 1])


%% sbesseljp2

x = linspace(0,20,1000)';
n = 0:5;
[X N] = ndgrid(x,n);
j1 = sbesseljp2(N,X);

figure(1),clf,hold all
plot(x,j1)
hold off
myplot('sbesseljp2(x)','x',' ')
ylim([-2 2])


%% sbesselh

x = linspace(0,20,1000)';
n = 0:5;
[X N] = ndgrid(x,n);
h1 = sbesselh(N,X);

figure(1),clf,hold all
plot(x,real(h1))
plot(x,imag(h1),'--')
hold off
myplot('sbesselh(x)','x',' ')
ylim([-2 2])


%% sbesselhp

x = linspace(0,20,1000)';
n = 0:5;
[X N] = ndgrid(x,n);
h1 = sbesselhp(N,X);

figure(1),clf,hold all
plot(x,real(h1))
plot(x,imag(h1),'--')
hold off
myplot('sbesselhp(x)','x',' ')
ylim([-2 2])

%% sbesselhp2

x = linspace(0,20,1000)';
n = 0:5;
[X N] = ndgrid(x,n);
h1 = sbesselhp2(N,X);

figure(1),clf,hold all
plot(x,real(h1))
plot(x,imag(h1),'--')
hold off
myplot('sbesselhp2(x)','x',' ')
ylim([-2 2])



