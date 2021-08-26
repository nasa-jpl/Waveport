function [ind] = indAzBz(l,lp,m,Lp)
% Row-major linear index of the sparse vector axial translation matrix
% 
% l,lp,m:   Harmonix index (l,lp,m)
% Lp:       Maximum column degree lp=1:Lp
%
% ind:      Sparse matrix linear index
%
% Dependencies: NAzBz

if m > min(l,lp) || l < 1 || lp < 1
    disp('bad index')
else
    if l == 1 
        ind1 = 0;
    else
        ind1 = NAzBz(l-1,Lp);
    end    
    if lp == 1
        ind2 = 0;        
    elseif l+1>=lp && ((l==1 && lp>=2) || (lp==2 && l>1))
        ind2 = 3;
    elseif l+1<lp && ((l==1 && lp>=2) || (lp==2 && l>1))
        ind2 = 2*l*(lp-2) + lp + 1;
    elseif l > 1 && (l + 1 >= lp) && l > 2
        ind2 = lp^2 - 1;        
    else
        ind2 = -l^2 + l + (2*l + 1)*(lp - 1); 
    end
    ind3 = min(l,lp) + m + 1;
    ind = ind1 + ind2 + ind3;
end
ind = round(ind);