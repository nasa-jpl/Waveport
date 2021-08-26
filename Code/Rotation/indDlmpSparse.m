function [ind] = indDlmpSparse(l,m,p,str)
% Returns linear index for sparse rotation matrix Dlmp
%
% l,m,p:    angular harmonic indices
% str:      'mono' for monopole
%
% ind:      sparse matrix linear index

if nargin == 4 && strcmp(str,'mono')
    if l < 0 || abs(m) > l
        error('bad index')
    end    
    if l == 0
        ind = 1;
        return
    else
        c1 = (2*l+1)*(l+p);
        c2 = l+m+1;
        ind = NDlmpSparse(l-1,str) + round(c1 + c2);
    end
else
    if l < 1 || abs(m) > l
        error('bad index')
    end
    c1 = (2*l+1)*(l+p);
    c2 = l+m+1;
    if l == 1      
        ind = round(c1 + c2);
    else
        ind = NDlmpSparse(l-1) + round(c1 + c2);
    end                
end
