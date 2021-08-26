function [clm dlm A B] = vectorTranslation(alm,blm,L,Lp,k,xji,rgstr,normstr)
% Vector translation matrix
%
% alm,blm:      Expansion coefficients in frame i
% L,Lp:         Maximum row/col degrees of the translation matrix
% xji:          3x1 Cartesian vector the points from frame i to frame j
% rgstr:        Use 'rg' for regular form of translation matrix
% normstr:      Use 'norm' for normalized vector wave functions
% 
% clm,dlm:      Expansion coefficients in frame j
% A,B:          [optional] Output the pair of full translation matrices
%
% Dependencies: cart2sph,AzBzSparse,DlmpSparseFast,NDlmpSparse

NL = L^2 + 2*L;
NLp = Lp^2 + 2*Lp;
Lmax = max([L Lp]); 
Lmin = min([L Lp]);
NLmax = max([NL NLp]);
NLmin = min([NL NLp]);
alm = alm(:);
blm = blm(:);
if length(alm) ~= NLp || length(blm) ~= NLp
    error('alm/blm array lengths do not match Lp')
end

% convert translation vector to spherical 
[rji thetaji phiji] = cart2sph(xji(1),xji(2),xji(3));

% compute sparse axial translation matrices
[row col Az Bz] = AzBzSparse(L,Lp,k,rji,rgstr,normstr);
Az = sparse(row,col,Az,NL,NLp);
Bz = sparse(row,col,Bz,NL,NLp);

% compute rotation matrix for larger of L and Lp
[row col D] = DlmpSparseFast(Lmax,phiji+pi/2,thetaji,0);
D1 = sparse(row,col,D,NLmax,NLmax);

% grab the matrix entries for the lesser of L and Lp
Ngrab = NDlmpSparse(Lmin);
D2 = sparse(row(1:Ngrab),col(1:Ngrab),D(1:Ngrab),NLmin,NLmin);

% mutliply depending on L and Lp
if Lp >= L
    tmpalm = D1*alm;
    tmpblm = D1*blm;
    clm = D2'*(Az*tmpalm + Bz*tmpblm);
    dlm = D2'*(Bz*tmpalm + Az*tmpblm);
else
    tmpalm = D2*alm;
    tmpblm = D2*blm;
    clm = D1'*(Az*tmpalm + Bz*tmpblm);
    dlm = D1'*(Bz*tmpalm + Az*tmpblm);
end

if nargout == 4
    if Lp >= L
        A = D2'*Az*D1;
        B = D2'*Bz*D1;
    else
        A = D1'*Az*D2;
        B = D1'*Bz*D2;
    end
end        