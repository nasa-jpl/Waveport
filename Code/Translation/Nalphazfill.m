function Nfill = Nalphazfill(L,Lp)
% Number of non-zero elements in the fill region of the 
% sparse scalar axial translation matrix.  Assumes L >= Lp
%
% L:    Maximum row order in actual matrix, l=0:L
% Lp:   Maximum column order in actual matrix lp=0:Lp
%
% Nfill:   Number of nonzero elements in the fill region

if L < Lp
    disp('bad index')
    return
end
Nfill = 1/6*Lp*(2*Lp^2 + 3*Lp + 1);
Nfill = round(Nfill);


