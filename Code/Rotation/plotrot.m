function plotrot(R,v,s,colorstr,linestr)
% Plot X', Y', Z' axes of a rotated coordinate system at location v.
%
% R:            3x3 rotation matrix
% v:            [x,y,z] location of coordinate frame
% s:            scale factor
% colorstr:     3x1 char array with axes colors
% linestr:      linestyle type

if nargin <= 3;
    colorstr = 'rbg';
    linestr = '-';
elseif nargin >= 4 && isempty(colorstr)
    colorstr = 'rbg';
elseif nargin >= 4 && length(colorstr) == 1;
    colorstr = repmat(colorstr,1,3);
end
if nargin == 5 && isempty(linestr)
    linestr = '-';
end
hold on
quiver3(v(1),v(2),v(3),R(1,1),R(2,1),R(3,1),s,'color',colorstr(1),'linestyle',linestr)
quiver3(v(1),v(2),v(3),R(1,2),R(2,2),R(3,2),s,'color',colorstr(2),'linestyle',linestr)
quiver3(v(1),v(2),v(3),R(1,3),R(2,3),R(3,3),s,'color',colorstr(3),'linestyle',linestr)
hold off
