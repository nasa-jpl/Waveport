function drawaxes(cen,lims)

Z = [0 0];
hold on
plot3([lims(1) lims(2)]+cen(1),Z+cen(2),Z+cen(3),'k')
plot3(Z+cen(1),[lims(3) lims(4)]+cen(2),Z+cen(3),'k')
plot3(Z+cen(1),Z+cen(2),[lims(5) lims(6)]+cen(3),'k')
hold off