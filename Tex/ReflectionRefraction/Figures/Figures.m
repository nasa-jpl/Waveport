
direc = '/Users/mshaynes/Desktop/Book/Waveport/Tex/ReflectionRefraction/Figures/';


%% Reflection on a plane

N = 5;
rng(7);
x1 = 10*rand(N,1);
y1 = 10*rand(N,1);
z1 = 10*rand(N,1);

x2 = 10*rand(N,1);
y2 = 10*rand(N,1);
z2 = 10*rand(N,1);

zo = 0;

[rx ry rz theta l1 l2] = reflectionPlane(x1,y1,z1,x2,y2,z2,zo);

h1 = figure(1),clf,hold all
scatter3(x1,y1,z1,'filled','b');
scatter3(x2,y2,z2,'filled','r');
scatter3(rx,ry,rz,'filled','g');
for n=1:N,
    plot3([x1(n) rx(n)]',[y1(n) ry(n)]',[z1(n) rz(n)]','b');
    plot3([x2(n) rx(n)]',[y2(n) ry(n)]',[z2(n) rz(n)]','b');
end
axis square
hold off
myplot('Reflection point on a plane','x','y','z')
view([1 1 0.6])


if 0
    saveimage(h1,[direc 'reflectionraysplane'],'epsc')
end



%% Reflection point on a sphere

r = 2;

nlat = 5;
[th phi] = discoball(nlat);
Npts = length(th);
rpts = 5;
[x2 y2 z2] = sph2cart(rpts*ones(size(th)),th,phi);

x1 = 4*ones(Npts,1);
y1 = -4*ones(Npts,1);
z1 = 4*ones(Npts,1);

[rx ry rz] = reflectionSphere(x1,y1,z1,x2,y2,z2,r);


[X Y Z] = sphere(40);
X = r*X;
Y = r*Y;
Z = r*Z;

h1 = figure(1),clf,hold all
scatter3(x1,y1,z1,'k','filled')
scatter3(x2,y2,z2,'filled')
scatter3(rx,ry,rz,'k','filled')
surf(X,Y,Z)
plot3([x1 rx]',[y1 ry]',[z1 rz]','b');
plot3([x2 rx]',[y2 ry]',[z2 rz],'r');
hold off
axis equal
shading flat
view([1 0.3 0.5])
myplot('Reflection Points on a Sphere','x','y','z',1)

if 0
    saveimage(h1,[direc 'reflectonraysphere'],'epsc')
end

%% Refraction at a plane


x1 = 0;
y1 = 0;
z1 = 5;
n1 = 1;
n2 = 2;
zo = 0;

N = 31;
x = linspace(-10,10,N);
y = 0;
[X2 Y2] = ndgrid(x,y);
Z2 = -5*ones(size(X2));


[RX RY RZ theta1 theta2 l1 l2] = refractionPlane(x1,y1,z1,X2,Y2,Z2,zo,n1,n2);

NN = numel(X2);
h1 = figure(1),clf,hold all
scatter3(x1,y1,z1,'filled')
scatter3(X2(:),Y2(:),Z2(:),'filled')
plot3([x1*ones(NN,1) RX(:)]',[y1*ones(NN,1) RY(:)]',[z1*ones(NN,1) RZ(:)]','b-');
plot3([X2(:) RX(:)]',[Y2(:) RY(:)]',[Z2(:) RZ(:)]','b-');
plot3([-12 12],[0 0],[0 0],'k')
hold off
view([0 -1 0])
myplot('Refraction Through a Plane','x','y','z')
zlim([-6 6])
xlim([-12 12])

if 0
    saveimage(h1,[direc 'refractionraysplane'],'epsc')
end



%% Refraction through circle



r = 1;
n1 = 1;
n2 = 2;

x1 = 0;
y1 = 2;
x2 = -0.5;
y2 = 0.1;

figure(1),clf,hold all
scatter(0,y1,'filled')
scatter(x2,y2,'filled')
plotcircle(0,0,r)
hold off
axis equal


[rx ry theta1 theta2 v1mag v2mag] = refractionCircle(x1,y1,x2,y2,n1,n2,r);

h1 = figure(1),clf,hold all
scatter(x1,y1,'filled')
scatter(x2,y2,'filled')
scatter([rx],[ry],'filled')
scatter(0,0,'filled')
plotcircle(0,0,r)
%plot([0 0],[0 y1])
plot([x2 rx],[y2 ry])
plot([x1 rx],[y1 ry])
%plot([0 rx],[0 ry])
hold off
axis equal
myplot('Refraction Into a Circle','x','y')
ylim([-1.1 2.1])

if 0
    saveimage(h1,[direc 'refractionrayscircle'],'epsc')
end

%% refraction though circle iterative


r = 1;
n1 = 1;
n2 = 2;

Npts = 20;
x2 = 0.1*linspace(-1,1,Npts)';
y2 = 0.8*ones(Npts,1);

x1 = 0*ones(Npts,1);
y1 = 1.2*ones(Npts,1);

figure(1),clf,hold all
scatter(x1,y1,'filled')
scatter(x2,y2,'filled')
plotcircle(0,0,r)
hold off
axis equal

nit = 40;
[rx ry theta1 theta2 thetal v1mag v2mag] = refractionCircleIterative(x1,y1,x2,y2,n1,n2,r);


h1=figure(1),clf,hold all
scatter(x1,y1,'filled')
scatter(x2,y2,'filled')
scatter([rx],[ry],'filled')
scatter(0,0,'filled')
plotcircle(0,0,r)
%plot([0 0],[0 y1])
plot([x2 rx]',[y2 ry]','b')
plot([x1 rx]',[y1 ry]','b')
%plot([0 rx],[0 ry])
hold off
axis equal
myplot('Refraction through a Circle','x','y')
xlim([-.2 .2])
ylim([0.75 1.25])


if 0
    saveimage(h1,[direc 'refractionrayscircleiterative'],'epsc')
end


%% refraction through a sphere

% sphere radius and media index of refractions
r = 1;
n1 = 1;
n2 = 2;

if 1 % points example 1
% define multiple interior points
    Npts = 5;
    x2 = 0.4*linspace(-1,1,Npts)';
    y2 = x2;
    z2 = 0.2*ones(Npts,1);
    [x2 y2 z2] = ndgrid(x2,y2,z2);
    
    % pair the same source point with all interior points
    sz = size(x2);
    x1 = 0*ones(sz);
    y1 = 0*ones(sz);
    z1 = 2*ones(sz);
else % points example 2
 % define multiple interior points
    Npts = 3;
    x2 = 0.1*linspace(-1,1,Npts)';
    y2 = x2;
    z2 = 0.8*ones(Npts,1);
    [x2 y2 z2] = ndgrid(x2,y2,z2);
    
    % pair the same source point with all interior points
    sz = size(x2);
    x1 = 0*ones(sz);
    y1 = 0*ones(sz);
    z1 = 1.2*ones(sz);
end


if 1
    type = 'dir';
    nit = [];
    typelab = 'Direct Solution';
else % iterative solution
    type = 'it';
    nit = [];
    typelab = 'Iterative Solution';
end

% compute sphere refraction
[rx ry rz theta1 theta2 thetal v1mag v2mag] = refractionSphere(x1,y1,z1,x2,y2,z2,n1,n2,r,type,nit);

% plot
h1 = figure(1),clf,hold all
scatter3(x1(:),y1(:),z1(:),'filled')
scatter3(x2(:),y2(:),z2(:),'filled')
scatter3(rx(:),ry(:),rz(:),'filled')
[X Y Z] = sphere(30);
surf(X,Y,Z,'FaceAlpha',0.4), shading flat
plot3([x2(:) rx(:)]',[y2(:) ry(:)]',[z2(:) rz(:)]','b')
plot3([x1(:) rx(:)]',[y1(:) ry(:)]',[z1(:) rz(:)]','b')
hold off
axis equal
myplot({'Refraction through a Sphere';typelab},'x','y','z')
view([0.3 1 0.5])


if 0
    saveimage(h1,[direc 'refractionspheredirect'],'epsc')
end

