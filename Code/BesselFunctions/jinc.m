function jc = jinc(x)
% jinc (Sombrero) function
%
% jc:   jinc function: jinc(x) = J_1(x)/x

jc = besselj(1,x)./x;
thresh = 1e-5;
jc((x<thresh)) = 1/2-x((x<thresh)).^2/16;