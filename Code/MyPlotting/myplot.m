function myplot(t,x,y,z,lwflag)

set(gcf,'color','w')
set(gca,'fontsize',14)
if nargin <= 4
    h2 = get(gca, 'Children');
    if ~isempty(h2)
        if isprop(h2,'linewidth')
            set(h2,'linewidth',2);
        end
    end
end
set(gca,'linewidth',1)
if nargin > 0
    title(t, 'interpreter', 'tex')
    xlabel(x);
    ylabel(y);
end
if nargin >= 4
   zlabel(z);
end
grid minor
if exist('dock.m') == 2
    dock
end
end

