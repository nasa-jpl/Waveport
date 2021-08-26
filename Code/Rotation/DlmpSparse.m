function [row col Dlmp] = DlmpSparse(L,alpha,beta,gamma,str)
% Rotation addition theorem for spherical harmonics, sparse matrix
%
% Y_l^m(theta,phi) = sum_{l,p} D_{l,m,p} Y_l^p(theta',phi')
%
% l = 1:L, for use with vector wave functions (no monopole).  
% (theta', phi') are angles seen from the rotated system.
% Rotation defined from unprimed to primed given by ZXZ 
% Euler angles alpha, beta, gamma 
%
% L:            Largest harmonic l, with l = 1:L, or l = 0:L
% alpha:        Z Euler angle (rad)
% beta:         X Euler angle (rad)
% gamma:        Z Euler angle (rad)
% str:          'mono' for monopole
%
% row,col:      row and column indices of matrix entries
% Dlmp:         Sparse rotation matrix entries
%
% Dependencies:     euler2rot,NDlmpSparse,indDlmpSparse,lm2ind,a_lmp,b_lmp,c_lmp,d_lmp

if L == 0
    Dlmp = 1;
    return
elseif L < 1
    disp('bad L')
    return
end
[Rot] = euler2rot(alpha,beta,gamma);
Rot = Rot.';
rxx = Rot(1,1);
rxy = Rot(1,2);
rxz = Rot(1,3);
ryx = Rot(2,1);
ryy = Rot(2,2);
ryz = Rot(2,3);
rzx = Rot(3,1);
rzy = Rot(3,2);
rzz = Rot(3,3);
F = [(ryy+rxx)/2 rxz/sqrt(2) (ryy-rxx)/2;
    rzx/sqrt(2) rzz -rzx/sqrt(2);
    (ryy-rxx)/2 -rxz/sqrt(2) (ryy+rxx)/2];
G = [(ryx-rxy)/2 ryz/sqrt(2) -(ryx+rxy)/2;
    -rzy/sqrt(2) 0 -rzy/sqrt(2);
    (ryx+rxy)/2 ryz/sqrt(2) (rxy-ryx)/2];
D1 = F+1i*G;

NL = NDlmpSparse(L);
Dlmp = zeros(NL,1);
row = zeros(NL,1);
col = zeros(NL,1);

l=1;
for m=-l:l,
for p=-l:l,
    inm = lm2ind(l,m); 
    inp = lm2ind(l,p); 
    inS = indDlmpSparse(l,m,p);
    row(inS) = inm;
    col(inS) = inp;
    Dlmp(inS) = D1(inm,inp);
end
end

for l=2:L,      
    for m = (-l+1):(l-1),
    for p = (-l+1):(l-1),    
        a1 = a_lmp(l,m,p);
        b1 = b_lmp(l,m,p);
        b2 = b_lmp(l,-m,p);        
        if a1~=0            
            tmp1 = a1*D1(2,2)*Dlmp(indDlmpSparse(l-1,m,p));
        else tmp1 = 0; end
        if b1~=0
            tmp2 = b1*D1(3,2)*Dlmp(indDlmpSparse(l-1,m-1,p));
        else tmp2 = 0; end
        if b2~=0
            tmp3 = b2*D1(1,2)*Dlmp(indDlmpSparse(l-1,m+1,p));
        else tmp3 = 0; end
        inS = indDlmpSparse(l,m,p);
        row(inS) = lm2ind(l,m); 
        col(inS) = lm2ind(l,p); 
        Dlmp(inS) = tmp1+tmp2+tmp3;
    end
    end
    m = -l;
    for p = (-l+1):(l-1),    
        b2 = b_lmp(l,-m,p);
        inS = indDlmpSparse(l,m,p);
        row(inS) = lm2ind(l,m);
        col(inS) = lm2ind(l,p);        
        Dlmp(inS) = b2*D1(1,2)*Dlmp(indDlmpSparse(l-1,m+1,p));        
    end
    m = l;
    for p = (-l+1):(l-1),
        b1 = b_lmp(l,m,p);
        inS = indDlmpSparse(l,m,p);
        row(inS) = lm2ind(l,m); 
        col(inS) = lm2ind(l,p);         
        Dlmp(inS) = b1*D1(3,2)*Dlmp(indDlmpSparse(l-1,m-1,p)); 
    end
    p = -l;
    for m = (-l):(l),
        c1 = c_lmp(l,m,-p);
        d1 = d_lmp(l,m,-p);
        d2 = d_lmp(l,-m,-p);
        if c1~=0
            tmp1 = c1*D1(2,1)*Dlmp(indDlmpSparse(l-1,m,p+1));
        else tmp1 = 0; end
        if d1~=0
            tmp2 = d1*D1(3,1)*Dlmp(indDlmpSparse(l-1,m-1,p+1));
        else tmp2 = 0; end
        if d2~=0
            tmp3 = d2*D1(1,1)*Dlmp(indDlmpSparse(l-1,m+1,p+1));
        else tmp3 = 0; end 
        inS = indDlmpSparse(l,m,p);
        row(inS) = lm2ind(l,m); 
        col(inS) = lm2ind(l,p); 
        Dlmp(inS) = tmp1+tmp2+tmp3;   
    end
    p = l;
    for m = (-l):(l),
        c1 = c_lmp(l,m,p);
        d1 = d_lmp(l,m,p);
        d2 = d_lmp(l,-m,p);       
        if c1~=0
             tmp1 = c1*D1(2,3)*Dlmp(indDlmpSparse(l-1,m,p-1));
        else tmp1 = 0; end
        if d1~=0
             tmp2 = d1*D1(3,3)*Dlmp(indDlmpSparse(l-1,m-1,p-1));
        else tmp2 = 0; end
        if d2~=0
            tmp3 = d2*D1(1,3)*Dlmp(indDlmpSparse(l-1,m+1,p-1));
        else tmp3 = 0; end

        inS = indDlmpSparse(l,m,p);
        row(inS) = lm2ind(l,m); 
        col(inS) = lm2ind(l,p); 
        Dlmp(inS) = tmp1+tmp2+tmp3;
    end 
end
if nargin == 5 && strcmp(str,'mono')
    NL = NDlmpSparse(L,'mono');
    Dlmp2 = zeros(NL,1);
    Dlmp2(2:end) = Dlmp;
    Dlmp = Dlmp2;
    Dlmp(1) = 1;
    row2 = zeros(NL,1);
    row2(2:end) = row + 1;
    row = row2;
    row(1) = 1;
    col2 = zeros(NL,1);
    col2(2:end) = col + 1;
    col = col2;
    col(1) = 1;       
end
end

% Helper functions
function a = a_lmp(l,m,p)
    a = sqrt((l+m)*(l-m)/((l+p)*(l-p)));
end
function b = b_lmp(l,m,p)
    b = sqrt((l+m)*(l+m-1)/(2*(l+p)*(l-p)));
end
function c = c_lmp(l,m,p)
    c = sqrt(2*(l+m)*(l-m)/((l+p)*(l+p-1)));
end
function d = d_lmp(l,m,p)
    d = sqrt((l+m)*(l+m-1)/((l+p)*(l+p-1)));
end

