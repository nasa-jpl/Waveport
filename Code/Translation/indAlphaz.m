function [ind] = indAlphaz(l,lp,m,Lp)
% Row-major linear index of the sparse scalar axial translation matrix
% 
% l,lp,m:   Harmonix index (l,lp,m)
% Lp:       Maximum column degree lp=0:Lp
%
% ind:      Sparse matrix linear index
%
% Dependencies: Nalaphz

if m > min(l,lp) || l < 0 || lp < 0
    disp('bad index')
else
    if l == 0 
        ind1 = 0;
    else
        ind1 = Nalphaz(l-1,Lp);
    end    
    if lp == 0
        ind2 = 0;
    elseif lp == 1 || ((l+1>lp) && (lp > 1))
        ind2 = lp^2;
    else
        ind2 = -l^2 - l + 2*l*lp + lp; 
    end
    ind3 = min(l,lp) + m + 1;
    ind = ind1 + ind2 + ind3;
end
ind = round(ind);


