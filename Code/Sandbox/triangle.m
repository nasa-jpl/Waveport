function [tr] = triangle(x,a,b,xo,yo)
% Generalized triangle function
%
% x:    range
% a:    half-width at base
% b:    vertical extent
% xo:   horizontal shift
% yo:   vertical shift
%
% tr:   
 
tr = max(b*(1-abs((x-xo)/a))+yo,yo);


