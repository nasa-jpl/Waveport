
direc = '/Users/mshaynes/Desktop/Work/Book/Chap6/';

set(h1,'PaperPositionMode','auto') 
saveas(h1,[direc 'TLinterp'],'epsc')

% 
% set(h2,'PaperPositionMode','auto') 
% saveas(h2,[direc 'rough1_corr'],'epsc')

set(h6,'PaperPositionMode','auto') 
saveas(h6,[direc 'samples'],'epsc')


set(h2,'PaperPositionMode','auto') 
saveas(h2,[direc 'filt1'],'epsc')
set(h3,'PaperPositionMode','auto') 
saveas(h3,[direc 'filt2'],'epsc')
set(h4,'PaperPositionMode','auto') 
saveas(h4,[direc 'filt3'],'epsc')
set(h5,'PaperPositionMode','auto') 
saveas(h5,[direc 'filt4'],'epsc')



set(h7,'PaperPositionMode','auto') 
saveas(h7,[direc 'digvsL'],'epsc')


set(h8,'PaperPositionMode','auto') 
saveas(h8,[direc 'storgvsL'],'epsc')


set(h9,'PaperPositionMode','auto') 
saveas(h9,[direc 'octree'],'epsc')


set(h9,'PaperPositionMode','auto') 
saveas(h9,[direc 'layerNN'],'epsc')

set(h9,'PaperPositionMode','auto') 
saveas(h9,[direc 'layerInt'],'epsc')

set(h9,'PaperPositionMode','auto') 
saveas(h9,[direc 'cone'],'epsc')



%% scattering matirx plots


ltop_lam = 0.75;
[li di kd L] = fmmL(ltop_lam,1,1e-3);

Level = fmmLevel(length(L),li,di,L);

er1 = 1;
er2 = 3;
a = 0.25;

nhat = unit([0 0 1]);

k1 = 2*pi*sqrt(er1);
k2 = 2*pi*sqrt(er2);


types = {'S_{11}','S_{12}','S_{21}','S_{22}'};

for type = 1:4,
    
    S = sMatrixDisk(type,Level.Theta(:),Level.Phi(:),Level.Theta(:),Level.Phi(:),k1,er1,k2,er2,nhat,a);

    cx = [0 0.2];
    h1 = figure(1),clf,
    subplot(221),imagesc(abs(S.TT)),colorbar,caxis(cx);
    myplot2('|S_{\theta,\theta}|','k_i (index)','k_s (index)')
    axis square

    subplot(222),imagesc(abs(S.TP)),colorbar,caxis(cx);
    myplot2('|S_{\theta,\phi}|','k_i (index)','k_s (index)')
    axis square

    subplot(223),imagesc(abs(S.PT)),colorbar,caxis(cx);
    myplot2('|S_{\phi,\theta}|','k_i (index)','k_s (index)')
    axis square

    subplot(224),imagesc(abs(S.PP)),colorbar,caxis(cx);
    myplot2('|S_{\phi,\phi}|','k_i (index)','k_s (index)')
    axis square

    axes('Units','Normal');
    h = title(types(type),'fontsize',20);
    set(gca,'visible','off')
    set(h,'visible','on')

    pause
    
%     direc = '/Users/mshaynes/Desktop/Work/Book/Chap6/';
% 
%     set(h1,'PaperPositionMode','auto') 
%     saveas(h1,[direc types{type} 'Disk2'],'epsc')

end



