function [Dlmp] = Dlmp(L,alpha,beta,gamma,str)
% Rotation addition theorem for spherical harmonics 
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
% str:          [optional] 'mono' for monopole
%
% Dlmp:         NxN Block diagonal rotation matrix
%               N = L^2 + 2*L
%               (l,m) index = l^2 + l + m
%               ...monopole 
%               N = L^2 + 2*L + 1
%               (l,m) index = l^2 + l + m + 1
%
% Dependencies:     euler2rot,lm2ind,a_lmp,b_lmp,c_lmp,d_lmp

if L == 0
    Dlmp = 1;
    return
elseif L < 1
    error('bad L')
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
tot = L^2 + 2*L;
Dlmp = zeros(tot);
Dlmp(1:3,1:3) = D1;
for l=2:L,      
    for m = (-l+1):(l-1),
    for p = (-l+1):(l-1),    
        a1 = a_lmp(l,m,p);
        b1 = b_lmp(l,m,p);
        b2 = b_lmp(l,-m,p);        
        if a1~=0            
            tmp1 = a1*D1(2,2)*Dlmp(lm2ind(l-1,m),lm2ind(l-1,p));
        else tmp1 = 0; end
        if b1~=0
            tmp2 = b1*D1(3,2)*Dlmp(lm2ind(l-1,m-1),lm2ind(l-1,p));
        else tmp2 = 0; end
        if b2~=0
            tmp3 = b2*D1(1,2)*Dlmp(lm2ind(l-1,m+1),lm2ind(l-1,p));
        else tmp3 = 0; end
        inm = lm2ind(l,m); 
        inp = lm2ind(l,p); 
        Dlmp(inm,inp) = tmp1+tmp2+tmp3;
    end
    end
    m = -l;
    for p = (-l+1):(l-1),    
        b2 = b_lmp(l,-m,p);
        Dlmp(lm2ind(l,m),lm2ind(l,p)) = ...
            b2*D1(1,2)*Dlmp(lm2ind(l-1,m+1),lm2ind(l-1,p));
    end
    m = l;
    for p = (-l+1):(l-1),
        b1 = b_lmp(l,m,p);
        Dlmp(lm2ind(l,m),lm2ind(l,p)) = ...
            b1*D1(3,2)*Dlmp(lm2ind(l-1,m-1),lm2ind(l-1,p));        
    end
    p = -l;
    for m = (-l):(l),
        c1 = c_lmp(l,m,-p);
        d1 = d_lmp(l,m,-p);
        d2 = d_lmp(l,-m,-p);
        if c1~=0
            tmp1 = c1*D1(2,1)*Dlmp(lm2ind(l-1,m),lm2ind(l-1,p+1));
        else tmp1 = 0; end
        if d1~=0
            tmp2 = d1*D1(3,1)*Dlmp(lm2ind(l-1,m-1),lm2ind(l-1,p+1));
        else tmp2 = 0; end
        if d2~=0
            tmp3 = d2*D1(1,1)*Dlmp(lm2ind(l-1,m+1),lm2ind(l-1,p+1));
        else tmp3 = 0; end
        Dlmp(lm2ind(l,m),lm2ind(l,p)) = tmp1+tmp2+tmp3;
    end
    p = l;
    for m = (-l):(l),
        c1 = c_lmp(l,m,p);
        d1 = d_lmp(l,m,p);
        d2 = d_lmp(l,-m,p);       
        if c1~=0
            tmp1 = c1*D1(2,3)*Dlmp(lm2ind(l-1,m),lm2ind(l-1,p-1));
        else tmp1 = 0; end
        if d1~=0
            tmp2 = d1*D1(3,3)*Dlmp(lm2ind(l-1,m-1),lm2ind(l-1,p-1));
        else tmp2 = 0; end
        if d2~=0
            tmp3 = d2*D1(1,3)*Dlmp(lm2ind(l-1,m+1),lm2ind(l-1,p-1));
        else tmp3 = 0; end
        Dlmp(lm2ind(l,m),lm2ind(l,p)) = tmp1+tmp2+tmp3;
    end 
end
if nargin == 5 && ~isempty(str)
    if ~strcmp(str,'mono')
        error('bad string')
    end
    tot = L^2+2*L+1;
    dlmp2 = zeros(tot);
    dlmp2(2:end,2:end) = Dlmp;
    Dlmp = dlmp2;
    Dlmp(1) = 1;
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

