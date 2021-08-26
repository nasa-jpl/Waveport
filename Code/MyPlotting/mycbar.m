function mycbar(t,cx)

hh = colorbar;
title(hh,t,'fontsize',14);
%title(hh,t,'fontsize',14,'interpreter','latex');
set(hh,'fontsize',14);
if nargin > 1
    caxis(cx)
end


