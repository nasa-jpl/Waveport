function [ind] = indDlmpSparseFast(l,m,p)
% Returns linear index for sparse rotation matrix Dlmp
%
% l,m,p:    angular harmonic indices
% str:      'mono' for monopole
%
% ind:      sparse matrix linear index

% c1 = (2*l+1)*(l+p);
% c2 = l+m+1;
% ind = NDlmpSparse(l-1) + round(c1 + c2);
% ind = round(4/3*(l-1)^3 + 4*(l-1)^2 + 11/3*(l-1) + (2*l+1)*(l+p) + l+m+1);

ind = round(4*l^3/3 + 2*l^2 + 5*l/3 + 2*l*p + m + p);
           

