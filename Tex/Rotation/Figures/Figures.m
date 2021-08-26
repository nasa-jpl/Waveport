



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
h1 =figure(1),clf

subplot(131)
plotrot(R0,origin,sz,'k',':')
plotrot(R1,origin,sz,[],'--')
view(vw)
axis equal
axis([-1 1 -1 1 -1 1])
myplot('Z Rotation','x','y','z')

subplot(132)
plotrot(R1,origin,sz,'k',':')
plotrot(R2,origin,sz,[],'--')
view(vw)
axis equal
axis([-1 1 -1 1 -1 1])
myplot('X'' Rotation','x','y','z')

subplot(133)
plotrot(R2,origin,sz,'k',':')
plotrot(R3,origin,sz,[],'--')
view(vw)
axis equal
axis([-1 1 -1 1 -1 1])
myplot('Z'''' Rotation','x','y','z')


if 0
    direc = './';
    saveimage(h1,[direc 'ZXZ'],'epsc');
end

%% Dlmp

L = 4;
alpha = 0;
beta = linspace(0,pi,9);
gamma = 0;
tot = L^2 + 2*L;

h1=figure(1),clf,
cx = [0 1];
for n = 1:length(beta),
    D = Dlmp(L,alpha,beta(n),gamma);
    subplot(3,3,n),imagesc(abs(D)),caxis(cx)
    myplot2(['\beta = ' num2str(180/pi*beta(n),4) ' ^o'],'linear (l,p)','linear (l,m)')
    axis square
end


if 0
    direc = './';
    saveimage(h1,[direc 'Dlmp'],'epsc');
end


%%





