



%% Point Source Expansion Coefficients 

Nka = 300;
ka = linspace(0,50,Nka);

L = 50;

alms = zeros(Nka,L+1);
for l=0:L,
alm = sqrt((2*l+1))*sbesselj(l,ka);
alms(:,l+1) = alm;
end

p = abs(alms).^2;
%p = p/p(1,1);

h1 = figure(1),clf
imagesc(ka,0:L,p'),colorbar
mycbar('$$|\widetilde{a}_{lm}|^2$$',[0 1])
set(gca,'ydir','normal')
hold all
plot(ka,ka,'g--','linewidth',2)
hold off
myplot('Spherical Harmonic Content of a Scalar Point Source','Radial Location of Point Target, ka','l')
grid off


grab = 150;
Lmax = ceil(1.1*ka(grab)*(1 + 1/ka(grab)));
h2 = figure(2),clf,hold all
stem(0:L,p(grab,:))
plot(Lmax*[1 1],[0 1.1*max(p(grab,:))],'k--')
hold off
myplot({'Expansion Coefficients of a Scalar Point Source';...
    ['ka = ' num2str(ka(grab),2)]},'l',' ')
legend('Coefficients','ceil(1.1ka(1+1/ka))')
ylabel('$$|\widetilde{a}_{lm}|^2$$','interpreter','latex')


if 0
    direc = './';
    saveimage(h1,[direc 'sphharmcont'],'epsc');
	saveimage(h2,[direc 'sphharmcontex'],'epsc');

end

