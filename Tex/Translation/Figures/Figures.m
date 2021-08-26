



%% figures for diagonalization 

k = 1;
rji = 100;
L = 6;
Lp = 6;

alp = alphaz(L,Lp,k,rji);

A = alp;
A(abs(alp)>0) = 1;
A = ~A;
colormap(gray)

tot1 = L^2+2*L+1;
tot2 = Lp^2+2*Lp+1;

h1 = figure(1),clf
imagesc(A)
hold on
% for l=0:(L-1);
%     plot([0 tot2+1],(lm2ind(l,l,'mono')+0.5)*[1 1],'k')
% end
% for l=0:(Lp-1);
%     plot((lm2ind(l,l,'mono')+0.5)*[1 1],[0 tot1+1],'k')
% end
%myplot('Scalar Axial Translation Matrix','Column Block','Row Block')
axis square
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
hold off
grid off

d = Dlmp(L,pi/8,pi/8,pi/8,'mono');

D = d;
D(abs(d)>0) = 1;
D = ~D;
colormap(gray)

h1 = figure(1),clf
imagesc(D)
hold on
% for l=0:(L-1);
%     plot([0 tot2+1],(lm2ind(l,l,'mono')+0.5)*[1 1],'k')
% end
% for l=0:(Lp-1);
%     plot((lm2ind(l,l,'mono')+0.5)*[1 1],[0 tot1+1],'k')
% end
%myplot('Scalar Axial Translation Matirx','Column Block','Row Block')
axis square
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
hold off
grid off




h1 = figure(1),clf
imagesc(zeros(size(A)))
hold on
% for l=0:(L-1);
%     plot([0 tot2+1],(lm2ind(l,l,'mono')+0.5)*[1 1],'k')
% end
% for l=0:(Lp-1);
%     plot((lm2ind(l,l,'mono')+0.5)*[1 1],[0 tot1+1],'k')
% end
%myplot('Scalar Axial Translation Matirx','Column Block','Row Block')
axis square
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
hold off
grid off



%% main translation matrix figure


k = 1;
rji = 100;
L = 8;
Lp = 8;

alp = alphaz(L,Lp,k,rji);

h1 = figure(1),clf
imagesc(db(alp)),colorbar
myplot('Scalar Axial Translation Matrix','(l'',m'') Linear Index','(l,m) Linear Index')
mycbar('(dB)',[-200 -20])
grid off


if 0
    direc = './';
    saveimage(h1,[direc 'trans1'],'epsc');
end


%% blocked out figure

k = 1;
rji = 100;
L = 8;
Lp = 4;

alp = alphaz(L,Lp,k,rji);

A = alp;
A(abs(alp)>0) = 1;
A = ~A;
colormap(gray)

tot1 = L^2+2*L+1;
tot2 = Lp^2+2*Lp+1;

h1 = figure(1),clf
imagesc(A)
hold on
for l=0:L;
    plot([0 tot2+1],(lm2ind(l,l,'mono')+0.5)*[1 1],'k')
end
for l=0:Lp;
    plot((lm2ind(l,l,'mono')+0.5)*[1 1],[0 tot1+1],'k')
end
plot(0.5*[1 1],[0 tot1+1],'k')


myplot('Scalar Axial Translation Matirx','Column Block','Row Block')
axis equal
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
hold off
grid off


%% full scalar translation matrix figure


k = 1;
rji = [12 5 15];
L = 8;
Lp = 8;

alp = alpha(L,Lp,k,rji);

h1 = figure(1),clf
imagesc(abs(alp)),colorbar
myplot('Scalar Translation Matrix','(l'',m'') Linear Index','(l,m) Linear Index')
mycbar('abs',[0 0.1])
grid off
axis square

if 0
    direc = './';
    saveimage(h1,[direc 'transfull'],'epsc');
end



%% full vector translation matrix figure


k = 1;
rji = [12 5 15];
L = 8;
Lp = 8;
rgstr = [];

for n=1:2,
    switch n
        case 1
            normstr = [];
            normtitle = 'Partially Normalized';
        case 2
            normstr = 'norm';
            normtitle = 'Fully Normalized';
    end

    [A B] = AB(L,Lp,k,rji,rgstr,normstr);

    h1 = figure(1),clf
    imagesc(abs(A)),colorbar
    myplot({'Vector Translation Matrix: A';normtitle},'(l'',m'') Linear Index','(l,m) Linear Index')
    mycbar('abs',[0 0.1])
    grid off
    axis square

    h2 = figure(2),clf
    imagesc(abs(B)),colorbar
    myplot({'Vector Translation Matrix: B';normtitle},'(l'',m'') Linear Index','(l,m) Linear Index')
    mycbar('abs',[0 0.1])
    grid off
    axis square


    if 0
        direc = './';
        saveimage(h1,[direc 'transfullA' num2str(n)],'epsc');
        saveimage(h2,[direc 'transfullB' num2str(n)],'epsc');
    end

end






