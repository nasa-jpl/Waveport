function [blm alp] = scalarTranslation(alm,L,Lp,k,xji,rgstr)
% Scalar translation matrix
%
% alm:      Expansion coefficients in frame i
% L,Lp:     Maximum row/col degrees of the translation matrix
% rji:      3x1 Cartesian vector the points from frame i to frame j
% rgstr:    Use 'rg' for regular form of translation matrix
% 
% blm:      Expansion coefficients in frame j
% alp:      [optional] entire translation matrix
%
% Dependencies: cart2sph,alphazSparse,DlmpSparse,NDlmpSparse

NL = L^2 + 2*L + 1;
NLp = Lp^2 + 2*Lp + 1;
Lmax = max([L Lp]);
Lmin = min([L Lp]);
NLmax = max([NL NLp]);
NLmin = min([NL NLp]);
alm = alm(:);
if length(alm) ~= NLp
    error('alm array length does not match Lp')
end

% convert translation vector to spherical 
[rji thetaji phiji] = cart2sph(xji(1),xji(2),xji(3));

% compute sparse axial translation
[row col alpz] = alphazSparse(L,Lp,k,rji,rgstr);
alpz = sparse(row,col,alpz,NL,NLp);

% compute rotation matrix for larger of L and Lp
[row col D] = DlmpSparse(Lmax,phiji+pi/2,thetaji,0,'mono');
D1 = sparse(row,col,D,NLmax,NLmax);

% grab the matrix entries for the lesser of L and Lp
Ngrab = NDlmpSparse(Lmin,'mono');
D2 = sparse(row(1:Ngrab),col(1:Ngrab),D(1:Ngrab),NLmin,NLmin);

% mutliply depending on L and Lp
if Lp >= L
    blm = D2'*(alpz*(D1*alm));
else
    blm = D1'*(alpz*(D2*alm));
end

if nargout == 2
    if Lp >= L
        alp = D2'*(alpz*(D1));
    else
        alp = D1'*(alpz*(D2));
    end
end