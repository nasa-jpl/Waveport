function Nl = NAzBz(L,Lp)
% Number of non-zero elements in sparse vector axial translation matrix
%
% L:    Maximum row degree l=1:L
% Lp:   Maximum column degree lp=1:Lp
%
% Nl:   Number of nonzero elements

if L == 1 && Lp == 1
    Nl = 3;
elseif Lp == 1 && L > 1
    Nl = 3*L;
elseif L == 1 && Lp > 1
    Nl = 3*Lp;
elseif L > 1 && Lp > 1 && L > Lp
    Nl = -1/3*Lp^3 + L*Lp^2 + 2*L*Lp + 1/3*Lp;
elseif L > 1 && Lp > 1 && Lp > L
	Nl = -1/3*L^3 + Lp*L^2 + 2*L*Lp + 1/3*L;
elseif L > 1 && Lp > 1 && Lp == L
    Nl = L^2 + 5*L + 2/3*Lp^3 + Lp^2 - 14/3*Lp;
else
    disp('bad index')
    return
end
Nl = round(Nl);