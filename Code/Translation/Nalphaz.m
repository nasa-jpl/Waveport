function Nl = Nalphaz(L,Lp)
% Number of non-zero elements in sparse scalar axial translation matrix
%
% L:    Maximum row order l=0:L
% Lp:   Maximum column order lp=0:Lp
%
% Nl:   Number of nonzero elements

if L == 0 && Lp == 0
    Nl = 1;
elseif L == 0 && (Lp>0)
    Nl = Lp + 1;
elseif Lp == 0 && (L>0)
    Nl = L + 1;
elseif (L>0) && (Lp>0) && (Lp<L) 
    Nl = -Lp^3/3 + Lp^2*L + Lp*(2*L + 4/3) + L + 1;
elseif (L>0) && (Lp>0) && (Lp>L)
    Nl = -L^3/3 + L^2*Lp + L*(2*Lp + 4/3) + Lp + 1;
elseif (L>0) && (Lp>0) && (Lp==L)
    Nl = 2/3*L^3 + 2*L^2 + 7/3*L + 1;
else
    disp('bad index')
    return
end
Nl = round(Nl);


