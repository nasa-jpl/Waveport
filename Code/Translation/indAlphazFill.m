function [ind] = indAlphazFill(l,lp,m,L,Lp)
% Row-major linear index for the sparse scalar axial translation matrix
% including the fill region.  Assumes that the fill region is added along 
% rows and L >= Lp.
% 
% l,lp,m:   Harmonic index (l,lp,m)
% L:        Maximum row degree l=0:L
% Lp:       Maximum column degree lp=0:Lp
%
% ind:      Sparse matrix linear index including the fill region
%
% Dependencies: Nalaphz, indAlphaz

if m > min(l,lp) || l < 0 || lp < 0 || Lp > L
    disp('bad index')
    return
elseif l <= L && lp <= Lp
    ind = indAlphaz(l,lp,m,Lp);
elseif l > L    
    ind = Nalphaz(L,Lp) + lp^2 + lp + m + 1;
    if l > L + 1 
        ind = ind + 1/6*(l-L-1)*(2*l^2-l*(4*L+6*Lp+7) ... 
            +2*L^2+L*(6*Lp+7)+6*(Lp+1)^2);
    end
else
    disp('bad index')
    return
end


