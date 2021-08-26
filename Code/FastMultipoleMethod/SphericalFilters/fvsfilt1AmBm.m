function [Am, Bm] = fvsfilt1AmBm(L,K,str)
% Preperatory function for fvsfilt1, which computes the matricies
% needed for the fast vector filter. Number of points are set by quadrature
% [I,J] = [2L+1,L+1], and [P,Q] = [2K+1,K+1]
%
% L:        Lower degree harmonic (to interpolate up from, or filter down to) 
% K:        Upper degree harmonic (to interpolate up to, or filter down from)
% str:      'interp' for interpolation from degree L to K
%           'filter' for filtering from degree K to L
%
% Am,Bm:    [JxQxI] matrix when interpolating from L to K
%           [QxJxI] matrix when filtering from K to L
%
% Dependencies: computeAmBm

if L > K, error('error, L > K'), end
I = 2*L+1;
J = L+1;
P = 2*K+1;
Q = K+1;
[muj wj] = legpts(J);
[muk wk] = legpts(Q);
if strcmp(str,'interp') % interp, call computeAmBm as coded
    [Am Bm] = computeAmBm(L,J,Q,I,muj,wj,muk);
elseif strcmp(str,'filter') % filter, call computeAmBm with j/k swapped
    [Am Bm] = computeAmBm(L,Q,J,I,muk,wk,muj);   
else
    error('error, bad string')
end
end
 
function [Am Bm] = computeAmBm(L,J,Q,I,muj,wj,muk)
% Computed matricies Am Bm for parameters passed from fvsfiltmatprep
    muj = muj(:);
    wj = wj(:);
    muk = muk(:);    
    aj = -sqrt(1-muj.^2);
    ak = -sqrt(1-muk.^2);
    bj = 1i./aj;
    bk = 1i./ak;
    Pj = Plm(L,muj); % includes monopole
    Ppj = Plmp(L,muj);
    Pk = Plm(L,muk);
    Ppk = Plmp(L,muk);

    % compute blocks
    M1m = zeros(J,Q,I);
    M2m = M1m;
    M3m = M1m;
    M4m = M1m;
    for m = -L:L,
        if m >= 0
            ind2 = m + 1;
        else
            ind2 = I + m + 1;
        end    
        if m == 0
            lowerlim = 1;
        else
            lowerlim = abs(m);
        end
        for l = lowerlim:L,
            ind = l^2 + l + m + 1; % monopole indexing
            const = 1/(l*(l+1));
            M1m(:,:,ind2) = M1m(:,:,ind2) + const*((Ppj(ind,:)).')*(Ppk(ind,:));
            M2m(:,:,ind2) = M2m(:,:,ind2) + const*((Pj(ind,:)).')*(Pk(ind,:));
            M3m(:,:,ind2) = M3m(:,:,ind2) + const*((Pj(ind,:)).')*(Ppk(ind,:));
            M4m(:,:,ind2) = M4m(:,:,ind2) + const*((Ppj(ind,:)).')*(Pk(ind,:));
        end            
    end 
    Am = zeros(J,Q,I);
    Bm = Am;
    for m = -L:L,
        if m >= 0
            ind2 = m + 1;
        else
            ind2 = I + m + 1;
        end   
        Am(:,:,ind2) = ((wj.*aj)*(ak.')).*M1m(:,:,ind2) - m^2*((wj.*bj)*(bk.')).*M2m(:,:,ind2);
        Bm(:,:,ind2) = m*((wj.*bj)*(ak.')).*M3m(:,:,ind2) + m*((wj.*aj)*(bk.')).*M4m(:,:,ind2);
    end 
end