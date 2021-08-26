

k = 2*pi;
xji = [10 10 40];
L = 3;
Lp = 4;
rgstr = 'rg';
tot = Lp^2 + 2*Lp + 1;
alm = ones(tot,1);

[blm alp] = scalarTranslation(alm,L,Lp,k,xji,rgstr);

figure(1),clf
imagesc(real(alp))

figure(2),clf
imagesc(imag(alp))


[a2] = alpha(L,Lp,k,xji,rgstr);

figure(3),clf
imagesc((abs(a2)>0))

figure(3),clf
imagesc(real(a2))

figure(4),clf
imagesc(imag(a2))


figure(1),clf
imagesc(real(alp))

figure(2),clf
imagesc(imag(alp))




%figure(5),hold all
%plot([real([alp(:) a2(:)]) imag([alp(:) a2(:)])],'o')

figure(5),clf,hold all
plot(real(alp(:)),'o')
plot(real(a2(:)),'.')
plot(imag(alp(:)),'s')
plot(imag(a2(:)),'x')
hold off




figure(6)
plot([real([alp(:)-a2(:)]) imag([alp(:)-a2(:)])])


%figure(7),clf
%imagesc(sign(real(alp))-sign(real(a2))),colorbar

%figure(8),clf
%imagesc(sign(imag(alp))-sign(imag(a2))),colorbar




