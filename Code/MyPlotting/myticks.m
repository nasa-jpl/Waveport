function myticks(xticks,yticks)
% custom function to quickly add tick lables on 2D plots

ticks = xticks;
tickslabels = [];
for m=1:length(ticks),tickslabels{m} = num2str(ticks(m));end
set(gca,'XTick',ticks)
set(gca,'XTickLabel',tickslabels)

ticks = yticks;
tickslabels = [];
for m=1:length(ticks),tickslabels{m} = num2str(ticks(m));end
set(gca,'YTick',ticks)
set(gca,'YTickLabel',tickslabels)