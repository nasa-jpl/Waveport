




%% momGmatrix2D

k = 2*pi;

% create a dense square grid of points
N = 20;
x = linspace(-2,2,N);
y = linspace(-2,2,N);
[X Y] = meshgrid(x,y);
dx = x(2)-x(1);

% return full matrix, set primed and unprimed coordinates to be the same
Xp = X;
Yp = Y;

% build the matrix
G = momGmatrix2D(dx,k,X,Y,Xp,Yp);

% plot
h1 = figure(1);,clf
imagesc(abs(G))
myplot('MoM Green''s Function Matrix - 2D','Point index (x'',y'')','Point index (x,y)')
colorbar
mycbar('abs(G)',[0 1e-2])
grid off

h2 = figure(2);,clf
imagesc(abs(G))
myplot('MoM Green''s Function Matrix - 2D','Point index (x'',y'')','Point index (x,y)')
colorbar
mycbar('abs(G)',[0 1e-2])
xlim([1 50])
ylim([1 50])
grid off

if 0
    direc = './';
    saveimage(h1,[direc 'green2d'],'epsc');
    saveimage(h2,[direc 'green2d_v2'],'epsc');
end


%% momGmatrix3D

k = 2*pi;

% create a dense square grid of points
N = 10;
x = linspace(-1,1,N);
y = linspace(-1,1,N);
z = linspace(-1,1,N);
[X Y Z] = meshgrid(x,y,z);
dx = x(2)-x(1);

% return full matrix, set primed and unprimed coordinates to be the same
Xp = X;
Yp = Y;
Zp = Z;

% build the matrix
G = momGmatrix3D(dx,k,X,Y,Z,Xp,Yp,Zp);

% plot
h1 = figure(1),clf
imagesc(abs(G))
myplot('MoM Green''s Function Matrix - 3D','Point index (x'',y'',z'')','Point index (x,y,z)')
colorbar
mycbar('abs(G)',5*[0 1e-2])
grid off

h2 = figure(2),clf
imagesc(abs(G))
myplot('MoM Green''s Function Matrix - 3D','Point index (x'',y'',z'')','Point index (x,y,z)')
colorbar
mycbar('abs(G)',5*[0 1e-2])
xlim([1 150])
ylim([1 150])
grid off

if 0
    direc = './';
    saveimage(h1,[direc 'green3d'],'epsc');
    saveimage(h2,[direc 'green3d_v2'],'epsc');
end



%% momGmatrixDyadic

k = 2*pi;

% create a dense square grid of points
N = 10;
x = linspace(-1,1,N);
y = linspace(-1,1,N);
z = linspace(-1,1,N);
[X Y Z] = meshgrid(x,y,z);
dx = x(2)-x(1);

% to return the full matrix, set primed and unprimed coordinates to be the same
Xp = X;
Yp = Y;
Zp = Z;

% build the matrix blocks
[Gxx Gyy Gzz Gxy Gxz Gyz] = momGmatrixDyadic(dx,k,X,Y,Z,Xp,Yp,Zp);

% plotting
for type = 1:6,
    switch type
        case 1
            tmp = Gxx; lab = 'Gxx'; cx = 3e-2*[0 1];
        case 2
            tmp = Gyy; lab = 'Gyy'; cx = 3e-2*[0 1];
        case 3
            tmp = Gzz; lab = 'Gzz'; cx = 3e-2*[0 1];
        case 4
            tmp = Gxy; lab = 'Gxy'; cx = 1e-2*[0 1];
        case 5
            tmp = Gxz; lab = 'Gxz'; cx = 1e-2*[0 1];
        case 6
            tmp = Gyz; lab = 'Gyz'; cx = 1e-2*[0 1];
    end

    figure(1),clf
    imagesc(abs(tmp))
    myplot(['MoM Dyadic Green''s Function Matrix Block ' lab],'Point index (x'',y'',z'')','Point index (x,y,z)')
    colorbar
    mycbar('abs(G)',cx)
    axis square
    grid off

    if 0
        direc = './';
        saveimage(h1,[direc 'green3dyadic' num2str(type)],'epsc');
    end

end

