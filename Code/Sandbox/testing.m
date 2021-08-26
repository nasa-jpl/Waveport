





Nlat = 20;
[th phi] = discoball(Nlat);
N = length(th);
[x y y] = sph2cart(ones(N,1),th,phi);
r2 = [x y y];

r1 = [0 0 2];
a = 0.5;

[i1 i2 l1 l2 delta_ind] = lineIntersectSphere2(r1,r2,a);

scatter3(x,y,y)


%%

N = 200;
theta = linspace(0,2*pi,N)';
y = cos(theta);
x = sin(theta);
R = 5;
r2 = R*[x y zeros(size(x))];

a = 2;
r1 = [0 1 0];

[i1 i2 l1 l2 delta_ind] = lineIntersectSphere2(r1,r2,a);

if 0
    R1 = repmat(r1,N,1);
    figure(1),clf,hold all
    plotcircle(0,0,a,'b')
    plot(r1(:,1),r1(:,2),'.')
    plot(r2(:,1),r2(:,2),'.')
    plot([i1(:,1) R1(:,1)]',[i1(:,2) R1(:,2)]','b')
    plot([i2(:,1) R1(:,1)]',[i2(:,2) R1(:,2)]','g')

    hold off
    axis equal
end

figure(2),clf
plot(l1)
ylim([0 4])

Z = r1(2);
U = sqrt(R^2 - 2*R*Z*cos(theta) + Z^2);
V = R*Z*cos(theta)-Z^2;
L = -V./U + sqrt((V./U).^2 - (Z^2-a^2));

figure(3),clf
plot(theta,L)
ylim([0 4])


