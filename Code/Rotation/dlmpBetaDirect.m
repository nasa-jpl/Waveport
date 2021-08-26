function dlmp = dlmpBetaDirect(L,beta,str)
% Direct computation of 'little-d' rotation matrix
%
% L:        maximum harmonic L
% beta:     Euler angle about X
% str:      'mono' to include monopole
%
% dlmp:     'little-d' rotation matrix
%
% Dependencies: lm2ind

N = L^2 + 2*L;
dlmp = zeros(N,N);
for l=1:L,
for m=-l:l,
for p=-l:l,
    row = lm2ind(l,m);
    col = lm2ind(l,p);
    f1 = factorial(l+p);
    f2 = factorial(l-p);
    f3 = factorial(l+m);
    f4 = factorial(l-m);
    const = sqrt(f1*f2/(f3*f4));
    tmp = 0;
    arg1 = l+m;
    arg2 = l-m;
    for u=0:(l+p),
       k1 = l+p-u;
       k2 = u;
       if k1 >= 0 && k1 <= arg1 && k2 >= 0 && k2 <= arg2
           t1=nchoosek(l+m,l+p-u);
           t2=nchoosek(l-m,u);
           signm = ((-1)^u)*((1i)^(m-p));
           c1 = (cos(beta/2)).^(2*l+p-m-2*u);
           s1 = (sin(beta/2)).^(m-p+2*u);        
           tmp = tmp + signm*t1*t2*c1*s1;
       end
    end
    dlmp(row,col) = const*tmp;
end
end    
end

if nargin == 3 && strcmp(str,'mono')
    tot = L^2+2*L+1;
    dlmp2 = zeros(tot);
    dlmp2(2:end,2:end) = dlmp;
    dlmp = dlmp2;
    dlmp(1) = 1;
end
