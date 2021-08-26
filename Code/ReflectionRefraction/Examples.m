

%% Reflection on a plane

% number of points
N = 5;
rng(7);

% sources
x1 = 10*rand(N,1);
y1 = 10*rand(N,1);
z1 = 10*rand(N,1);

% observers
x2 = 10*rand(N,1);
y2 = 10*rand(N,1);
z2 = 10*rand(N,1);

% plane level
zo = 0;

% compute reflection points
[rx ry rz theta l1 l2] = reflectionPlane(x1,y1,z1,x2,y2,z2,zo);

% plot
figure(1),clf,hold all
scatter3(x1,y1,z1,'filled','b');
scatter3(x2,y2,z2,'filled','r');
scatter3(rx,ry,rz,'filled','g');
for n=1:N,
    plot3([x1(n) rx(n)]',[y1(n) ry(n)]',[z1(n) rz(n)]','b');
    plot3([x2(n) rx(n)]',[y2(n) ry(n)]',[z2(n) rz(n)]','b');
end
axis equal
hold off
myplot('Reflection point on a plane','x','y','z')
view([1 1 1])



%% Reflection point on a sphere

% sphere radius
r = 2;

% define observation point around the sphere
nlat = 10;
[th phi] = discoball(nlat);
Npts = length(th);
rpts = 3;
[x2 y2 z2] = sph2cart(rpts*ones(size(th)),th,phi);

% define same source point to match observation points
x1 = 4*ones(Npts,1);
y1 = -4*ones(Npts,1);
z1 = 4*ones(Npts,1);

% compute the reflection points
[rx ry rz] = reflectionSphere(x1,y1,z1,x2,y2,z2,r);

% create sphere boundary for plotting
[X Y Z] = sphere(40);
X = r*X;
Y = r*Y;
Z = r*Z;

% plot
figure(1),clf,hold all
scatter3(x1,y1,z1,'k','filled')
scatter3(x2,y2,z2,'filled')
scatter3(rx,ry,rz,'k','filled')
surf(X,Y,Z)
plot3([x1 rx]',[y1 ry]',[z1 rz]','b');
plot3([x2 rx]',[y2 ry]',[z2 rz],'r');
hold off
axis equal
shading flat
view([1 1 1])




%% Refraction at a plane


% single source point
x1 = 0;
y1 = 0;
z1 = 5;

% index of refrations for upper and lower media
n1 = 1;
n2 = 2;

% interface level
zo = 0;

example = 3;
switch example
    case 1  % line of subsurface points
        N = 31;
        x = linspace(-5,5,N);
        y = 0;
        [X2 Y2] = ndgrid(x,y);
        Z2 = -4*ones(size(X2));
    case 2 % grid of subsurface points include point above
        N = 7;
        x = linspace(-5,5,N);
        z = linspace(-5,1,N);
        [Z2 X2] = ndgrid(z,x);
        Y2 = 0*ones(size(X2));
    case 3 % horizontal plane of subsurface points
        N = 11;
        x = linspace(-5,5,N);
        y = x;
        [X2 Y2] = ndgrid(x,y);
        Z2 = -4*ones(size(X2));
end
    
% compute direct solution and time it
tic
[RX RY RZ theta1 theta2 l1 l2] = refractionPlane(x1,y1,z1,X2,Y2,Z2,zo,n1,n2);
toc

% compute iterative solution and time it
tic
[RX RY RZ theta1 theta2 l1 l2] = refractionPlane(x1,y1,z1,X2,Y2,Z2,zo,n1,n2,'it');
toc

% plote
NN = numel(X2);
figure(1),clf,hold all
scatter3(x1,y1,z1,'filled')
scatter3(X2(:),Y2(:),Z2(:),'filled')
plot3([x1*ones(NN,1) RX(:)]',[y1*ones(NN,1) RY(:)]',[z1*ones(NN,1) RZ(:)]','b-');
plot3([X2(:) RX(:)]',[Y2(:) RY(:)]',[Z2(:) RZ(:)]','b-');
plot3([0 0],[-5 5],[0 0],'k')
plot3([-5 5],[0 0],[0 0],'k')
hold off
view([1 1 0.3])
myplot()


%% Refraction through a circle

% sphere radius and media index of refractions
r = 1;
n1 = 1;
n2 = 2;

% exterior point
x1 = 1;
y1 = 1.5;

% interior point
x2 = 0.5;
y2 = 0;


% compute refraction
[rx ry theta1 theta2 thetal v1mag v2mag] = refractionCircle(x1,y1,x2,y2,n1,n2,r);

% plot
figure(1),clf,hold all
scatter(x1,y1,'filled')
scatter(x2,y2,'filled')
scatter([rx],[ry],'filled')
scatter(0,0,'filled')
plotcircle(0,0,r)
plot([x2 rx],[y2 ry])
plot([x1 rx],[y1 ry])
hold off
axis equal
myplot('Refraction through circle','x','y')


%% Refraction through a circle - iterative

% sphere radius and media index of refractions
r = 1;
n1 = 1;
n2 = 2;

% define multiple interior points
Npts = 21;
x2 = 0.1*linspace(-1,1,Npts)';
y2 = 0.8*ones(Npts,1);

% pair one source point with all interior points
x1 = 0*ones(Npts,1);
y1 = 1.3*ones(Npts,1);

% plot geometry
figure(1),clf,hold all
scatter(x1,y1,'filled')
scatter(x2,y2,'filled')
plotcircle(0,0,r)
hold off
axis equal


% compute iterative solution using default number of iterations
[rx ry theta1 theta2 thetal v1mag v2mag] = refractionCircleIterative(x1,y1,x2,y2,n1,n2,r);

% plot
figure(1),clf,hold all
scatter(x1,y1,'filled')
scatter(x2,y2,'filled')
scatter(rx,ry,'filled')
scatter(0,0,'filled')
plotcircle(0,0,r)
plot([x2 rx]',[y2 ry]','b')
plot([x1 rx]',[y1 ry]','b')
hold off
axis equal
myplot('Refraction through a Circle','x','y')
xlim([-.2 .2])
ylim([0.75 1.25])


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
figure(1),clf,hold all
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



