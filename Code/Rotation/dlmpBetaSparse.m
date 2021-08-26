function [row col dlmp] = dlmpBetaSparse(L,beta,str)
% Recursive algorithm for ZXZ 'little-d' complex rotation matrix
% computed on sparse matrix
%
% L:        maximum harmonic L
% beta:     Euler angle about X
% str:      'mono' to include monopole
%
% row,col:      row and column indices of matrix entries
% dlmp:         Sparse 'little-d' rotation matrix
%
% Dependencies: lm2ind,NDlmpSparse,indDlmpSparse,

% includes monopole for initial condition
NL = NDlmpSparse(L,'mono');
dlmp = zeros(NL,1);
row = zeros(NL,1);
col = zeros(NL,1);
dlmp(1,1) = 1;
x = beta/2;
for l=1:L,
    m=l;
    for p=-l:l
        A = Almp(l,m,p);
        B = Blmp(l,m,p);
        C = Clmp(l,m,p);        
        if A~=0            
            tmp1 = cos(x)^2*A*...
                dlmp(indDlmpSparse(l-1,m-1,p-1,'mono'));
        else tmp1 = 0; end
        if B~=0
            tmp2 = -2*sin(x)*cos(x)*B*...
                dlmp(indDlmpSparse(l-1,m-1,p,'mono'));
        else tmp2 = 0; end
        if C~=0
            tmp3 = sin(x)^2*C*...
                dlmp(indDlmpSparse(l-1,m-1,p+1,'mono'));
        else tmp3 = 0; end
        inS = indDlmpSparse(l,m,p,'mono');
        row(inS) = lm2ind(l,m,'mono'); 
        col(inS) = lm2ind(l,p,'mono'); 
        dlmp(inS) = (tmp1+tmp2+tmp3);
    end
    
    m=-l;
    for p=-l:l
        D = Dlmp(l,m,p);
        E = Elmp(l,m,p);
        F = Flmp(l,m,p);
        if D~=0            
            tmp1 = sin(x)^2*D*...
                dlmp(indDlmpSparse(l-1,m+1,p-1,'mono'));
        else tmp1 = 0; end
        if E~=0
            tmp2 = 2*sin(x)*cos(x)*E*...
                dlmp(indDlmpSparse(l-1,m+1,p,'mono'));
        else tmp2 = 0; end
        if F~=0
            tmp3 = cos(x)^2*F*...
                dlmp(indDlmpSparse(l-1,m+1,p+1,'mono'));
        else tmp3 = 0; end
        inS = indDlmpSparse(l,m,p,'mono');
        row(inS) = lm2ind(l,m,'mono'); 
        col(inS) = lm2ind(l,p,'mono'); 
        dlmp(inS) = (tmp1+tmp2+tmp3);
    end
    
    for m=(-l+1):(l-1),
    for p=-l:l
        G = Glmp(l,m,p);
        H = Hlmp(l,m,p);
        I = Ilmp(l,m,p);
        if G~=0            
            tmp1 = sin(x)*cos(x)*G*...
                dlmp(indDlmpSparse(l-1,m,p-1,'mono'));
        else tmp1 = 0; end
        if H~=0
            tmp2 = (cos(x)^2 - sin(x)^2)*H*...
                dlmp(indDlmpSparse(l-1,m,p,'mono'));
        else tmp2 = 0; end
        if I~=0
            tmp3 = -sin(x)*cos(x)*I*...
                dlmp(indDlmpSparse(l-1,m,p+1,'mono'));
        else tmp3 = 0; end
        inS = indDlmpSparse(l,m,p,'mono');
        row(inS) = lm2ind(l,m,'mono'); 
        col(inS) = lm2ind(l,p,'mono'); 
        dlmp(inS) = (tmp1+tmp2+tmp3);
    end
    end
        
    
end
       
% convert to complex
for l=1:L,
for m=-l:l,
for p=-l:l,
    inS = indDlmpSparse(l,m,p,'mono');
    dlmp(inS) = ((1i)^(m-p)*(-1)^(m-p))*dlmp(inS);    
end
end
end

% no monopole
if nargin == 2
   row = row(2:end)-1;
   col = col(2:end)-1;
   dlmp = dlmp(2:end);
end
end

function a = Almp(l,m,p)
    a = sqrt((l+p)*(l+p-1)/((l+m)*(l+m-1)));
end
function a = Blmp(l,m,p)
    a = sqrt((l+p)*(l-p)/((l+m)*(l+m-1)));
end
function a = Clmp(l,m,p)
    a = sqrt((l-p)*(l-p-1)/((l+m)*(l+m-1)));
end
function a = Dlmp(l,m,p)
    a = sqrt((l+p)*(l+p-1)/((l-m)*(l-m-1)));
end
function a = Elmp(l,m,p)
    a = sqrt((l+p)*(l-p)/((l-m)*(l-m-1)));
end
function a = Flmp(l,m,p)
    a = sqrt((l-p)*(l-p-1)/((l-m)*(l-m-1)));
end
function a = Glmp(l,m,p)
    a = sqrt((l+p)*(l+p-1)/((l-m)*(l+m)));
end
function a = Hlmp(l,m,p)
    a = sqrt((l+p)*(l-p)/((l-m)*(l+m)));
end
function a = Ilmp(l,m,p)
    a = sqrt((l-p)*(l-p-1)/((l-m)*(l+m)));
end



