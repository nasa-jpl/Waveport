function dlmp = dlmpBetaFast(L,beta,str)
% Faster recursive algorithm for ZXZ 'little-d' complex rotation matrix
%
% L:        maximum harmonic L
% beta:     Euler angle about X
% str:      'mono' to include monopole
%
% dlmp:     'little-d' rotation matrix
%
% Dependencies:

% includes monopole for initial condition
N = L^2 + 2*L + 1;
dlmp = zeros(N,N);
dlmp(1,1) = 1;
x = beta/2;
CSvec = [cos(x)^2; -2*sin(x)*cos(x); sin(x)^2;...
    sin(x)^2; 2*sin(x)*cos(x); cos(x)^2;...
    sin(x)*cos(x); (cos(x)^2 - sin(x)^2); -sin(x)*cos(x)];

for l=1:L,
    m=l;
    for p=-l:l
        A = sqrt((l+p)*(l+p-1)/((l+m)*(l+m-1)));
        B = sqrt((l+p)*(l-p)/((l+m)*(l+m-1)));
        C = sqrt((l-p)*(l-p-1)/((l+m)*(l+m-1)));
        inlm1mm1 = (l-1)^2 + l - 1 + m;
        inlm1p = (l-1)^2 + l + p;
        if A~=0            
            tmp1 = CSvec(1)*A*dlmp(inlm1mm1,inlm1p-1);
        else tmp1 = 0; end
        if B~=0
            tmp2 =  CSvec(2)*B*dlmp(inlm1mm1,inlm1p);
        else tmp2 = 0; end
        if C~=0
            tmp3 =  CSvec(3)*C*dlmp(inlm1mm1,inlm1p+1);
        else tmp3 = 0; end
        inm = l^2 + l + m + 1; 
        inp = l^2 + l + p + 1; 
        dlmp(inm,inp) = (tmp1+tmp2+tmp3);
    end
    
    m=-l;
    for p=-l:l
        D = sqrt((l+p)*(l+p-1)/((l-m)*(l-m-1)));
        E = sqrt((l+p)*(l-p)/((l-m)*(l-m-1)));
        F = sqrt((l-p)*(l-p-1)/((l-m)*(l-m-1)));
        inlm1mp1 = (l-1)^2 + l + m + 1;
        inlm1p = (l-1)^2 + l + p;
        if D~=0            
            tmp1 =  CSvec(4)*D*dlmp(inlm1mp1,inlm1p-1);
        else tmp1 = 0; end
        if E~=0
            tmp2 =  CSvec(5)*E*dlmp(inlm1mp1,inlm1p);
        else tmp2 = 0; end
        if F~=0
            tmp3 =  CSvec(6)*F*dlmp(inlm1mp1,inlm1p+1);
        else tmp3 = 0; end
        inm = l^2 + l + m + 1; 
        inp = l^2 + l + p + 1; 
        dlmp(inm,inp) = (tmp1+tmp2+tmp3);
    end
    
    for m=(-l+1):(l-1),
    for p=-l:l
        inlm1m = (l-1)^2 + l + m;
        inlm1p = (l-1)^2 + l + p;
        G = sqrt((l+p)*(l+p-1)/((l-m)*(l+m)));
        H = sqrt((l+p)*(l-p)/((l-m)*(l+m)));
        I = sqrt((l-p)*(l-p-1)/((l-m)*(l+m)));
        if G~=0            
            tmp1 =  CSvec(7)*G*dlmp(inlm1m,inlm1p-1);
        else tmp1 = 0; end
        if H~=0
            tmp2 =  CSvec(8)*H*dlmp(inlm1m,inlm1p);
        else tmp2 = 0; end
        if I~=0
            tmp3 =  CSvec(9)*I*dlmp(inlm1m,inlm1p+1);
        else tmp3 = 0; end
        inm = l^2 + l + m + 1; 
        inp = l^2 + l + p + 1; 
        dlmp(inm,inp) = (tmp1+tmp2+tmp3);
    end
    end
end

% convert to complex
for l=1:L,
    m=-l:l;
    p=-l:l;
    inm = l^2 + l + m + 1; 
    inp = l^2 + l + p + 1; 
    [M, P] = ndgrid(m,p);
    D = M-P;
    dlmp(inm,inp) = ((-1i).^(D)).*dlmp(inm,inp);    
end

% no monopole
if nargin == 2
   dlmp = dlmp(2:end,2:end);
end
end
