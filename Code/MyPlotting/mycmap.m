function mycmap(minind)

%map1 = colormap(gray);
map1 = colormap(parula);
%map1 = flipud(map1);
%minind = 10;
map2 = zeros(size(map1));
xi = linspace(minind,64,64);
for n=1:3,
   map2(:,n) = interp1((minind:64)',map1(minind:64,n),xi,'linear'); 
end
colormap(map2)