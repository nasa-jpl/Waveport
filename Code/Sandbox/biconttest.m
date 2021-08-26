
%% bicont


fv = 0.3;
lave = 5e-3;
kave = 1500;
lave = 2*pi/kave;
b = 10;

scale = kave/(b+1);
shape = b+1;
alp = erfinv(1 - 2*fv);

Nx = 251;
x = linspace(-0.05,0.05,Nx);
y = x;
z = 0;
[X Y Z] = meshgrid(x,y,z);

rng(1);
[Theta S] = bicont(X,Y,Z,lave,b,fv);

mean(S(:))
var(S(:))

figure(1),clf,
imagesc(x*1000,y*1000,S),colorbar
colormap(gray)
myplot('Fluctuating Field','x (mm)','y (mm)')


figure(1),clf,
imagesc(x*1000,y*1000,Theta)
colormap(gray)
myplot('Indicator Function','x (mm)','y (mm)')

mean(Theta(:))








%%


dx = x(2)-x(1);
F = fft2(Theta);
C = 1/numel(Theta)*real((ifft2(F.*conj(F))));
%C = conv2(T,T,'same');
%C = C./max(C(:));

figure(2),clf,hold all
plot((C(1,:)))

CC = fftshift(C(1,:));

r = abs(x)';

r = linspace(0,1e-3,1000)';
[G Calp Cs] = bicontcorr(r,lave,b,fv);


figure(3),clf,hold all
plot(r,[Cs Calp]) 
ylim([-0.2 1])

figure(4),clf,hold all
plot(abs(x),[CC(:) G(:)])
xlim([0 0.02])


%%




r = linspace(0,0.01,1000)';

kc = 1500;
b = 10;

phi = atan(kc*r/(b+1));

Cs = cos(phi).^(b+1).*(phi./sin(phi)).*(sin(b*phi)./(b*phi));

plot(r,Cs)
grid on

%%
N = 10^4;

fv = 0.3;
kave = 10;
b = 14;

scale = kave/(b+1);
shape = b+1;

k = linspace(0,50,1000)';
p = 1/gamma(b+1)*(b+1)^(b+1)/kave*(k/kave).^b.*exp(-(b+1)*k/kave);

r = gamrnd(shape,scale,N,1);

figure(1),clf,hold all
histogram(r,'normalization','pdf')
plot(k,p)
hold off


%% levin's t-transform

K = 7;
tnk = zeros(K+2,K+1);

p = 1:9;
an = -(-2).^p./p;
Sn = cumsum(an);

tnk(:,1) = Sn;

for k = 1:K,
    for n=1:(K-k+1),
        num = 0;
        den = 0;
        for j=0:k,
            coef = (-1)^j*nchoosek(k,j)*((n+j)/(n+k)).^(k-1);
            num = num + coef*Sn(n+j)/(Sn(n+j+1)-Sn(n+j));
            den = den + coef/(Sn(n+j+1)-Sn(n+j));
        end
        tnk(n,k+1) = num/den;
    end
end


for k = 1:K,
    for n=1:(K-k+1),
        tnk(n,k+1) = levint(Sn,k,n);
    end
end

levint(Sn,K,1)



