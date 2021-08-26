function h = plotcircle(x,y,r,colorstr)

if nargin<4
    colorstr = 'k'
end
th = 0:pi/1000:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit,colorstr,'linewidth',1);
