function [A B] = AB(L1,L2,k,rs)
% 3D vector addition theorem

% A and B expand outgoing spherical vector wave function from the origin
%  as incoming waves about rs.  Implementation of Chew's recursive algorithm 
% not dependent on the scalar translation matrix alpha.

% Dependencies: SBESSELH.m, subfunctions

% L1, L2: maximum order l in row and col
% k: background wave number
% rs: [x, y, z] vector from source frame to receiving frame, units (m)
% alpha: N1xN2 square translation matrix
%                N1 = L1^2 + 2*L1 + 1
%                N2 = L2^2 + 2*L2 + 1
%                orders l1 = 0:L1, and l2 = 0:L2, all m
%                (l,m) = l^2 + l + m + 1

% W. C. Chew and Y. M. Wang, "Efficient ways to compute 
% the vector addition theorem," J. Electromagn. Waves Applicat., 
% vol. 7, no. 5, pp. 651-665, 1993.

% Mark Haynes
% University of Michigan
% 2006-2011


i = complex(0,1);

% convert rs to spherical coordinates
x = rs(1);
y = rs(2);
z = rs(3);
r = sqrt(x.^2 + y.^2 + z.^2);
ct = cos(atan2(sqrt(x.^2+y.^2),z)); 
phi = atan2(y,x);

% find greater of L1 and L2
L = max([L1 L2]);

tot1 = L1^2 + 2*L1;     % final row dimension
tot2 = L2^2 + 2*L2;     % final col dimension
tot3 = (2*L-1)^2 + 2*(2*L-1);   % temp row dimension

% generate (l,m) index table 
tota = (2*L)^2 + 2*(2*L) + 1;
tab = zeros(tota,2);
cnt = 1;
for l=0:2*L,
for m=(-l):l,
    tab(cnt,1) = l;
    tab(cnt,2) = m;
    cnt = cnt+1;
end
end

A = zeros(tot3,tot2);
B = zeros(tot3,tot2);


% need alpha_vuoo up thougth 2*L+1, because we 
% must fill out the initial conditions through 2*L
% which require alpha_(v+1)u00
alpha_vuoo = zeros(tota,1);

% alpha_{vu,00}
for in1 = 1:tota,
    
    lp = tab(in1,1);
    mp = tab(in1,2);
     
    ex = exp(i*(-mp)*phi);
    leg = legendre(lp,ct);
    if abs(-mp) > lp
        leg = 0;
    elseif sign(-mp) == -1
        leg = ((-1)^(-mp))*(factorial(lp-abs(-mp))/...
            factorial(lp+abs(-mp)))*leg(abs(-mp)+1);
    else
        leg = leg(-mp+1);
    end

    h = sbesselh(lp,k*r);
    nl = sqrt((2*lp+1)*factorial(lp+mp)/4/pi/factorial(lp-mp));
    clm = sqrt(4*pi)*((-1)^(lp+mp)); 

    alpha_vuoo(in1,1) = clm*nl*leg*ex*h;
end

% initial conditions
% A_vu10, A_vu1pm1, B_vu10, B_vu1pm1
for v=1:(2*L-1),
for u=-v:v,           
    vp = v+1;
    vm = v-1;
    up = u+1;
    um = u-1;    
    vu = v^2 + v + u + 1;
    vu1 = vm^2 + vm + u + 1;
    vu2 = vp^2 + vp + u + 1;
    vu3 = vp^2 + vp + um + 1;
    vu4 = vp^2 + vp + up + 1;
    vu5 = vm^2 + vm + um + 1;
    vu6 = vm^2 + vm + up + 1;
    vu7 = v^2 + v + um + 1;
    vu8 = v^2 + v + up + 1;
        
    if vu1 < 1 || vm < 0 || u < -vm || u > vm, tmp1 = 0;
    else tmp1 = zp(vm,u)*alpha_vuoo(vu1); end
    
    if vu2 > tota || u < -vp || u > vp, tmp2 = 0;
    else tmp2 = zm(vp,u)*alpha_vuoo(vu2); end
    
    if vu3 > tota || um < -vp || um > vp, tmp3 = 0;
    else tmp3 = emp(vp,um)*alpha_vuoo(vu3); end
    
    if vu4 > tota || up < -vp || up > vp, tmp4 = 0;
    else tmp4 = emm(vp,up)*alpha_vuoo(vu4); end
    
    if vu5 < 1 || vm < 0 || um < -vm || um > vm, tmp5 = 0;
    else tmp5 = epp(vm,um)*alpha_vuoo(vu5); end
    
    if vu6 < 1 || vm < 0 || up < -vm || up > vm, tmp6 = 0;
    else tmp6 = epm(vm,up)*alpha_vuoo(vu6); end
    
    if vu7 < 1 || um < -v || um > v, tmp7 = 0; 
    else tmp7 = eop(v,um)*alpha_vuoo(vu7); end
    
    if vu8 > tota || up < -v || up > v, tmp8 = 0;
    else tmp8 = eom(v,up)*alpha_vuoo(vu8); end
          
    A(vu-1,2) = sqrt(3)*(tmp1 + tmp2);
    B(vu-1,2) = sqrt(3)*zo(v,u)*alpha_vuoo(vu);
    A(vu-1,3) = -sqrt(3/2)*(tmp3 + tmp5);
    A(vu-1,1) = sqrt(3/2)*(tmp4 + tmp6);
    B(vu-1,3) = -sqrt(3/2)*tmp7;    
    B(vu-1,1) = sqrt(3/2)*tmp8;

end
end

% n+1,n+1
for n=1:(L2-1),
    nn = n^2 + 2*n;
    np = n+1;
    npnp = np^2 + 2*np;
    gppnn = gpp(n,n);

    for v = 1:(2*L-n-1),
    for u = (-v):v,
        vp = v+1;
        vm = v-1;
        um = u-1;
        vu = v^2 + v + u;
        vu1 = vp^2 + vp + um;
        vu2 = vm^2 + vm + um;
        vu3 = v^2 + v + um;
        gvpum = gmp(vp,um);
        gvmum = gpp(vm,um);
        gvum = gop(v,um);       

        if vu1 > tot3  || vu1 < 1 || um < -vp || um > vp
            tmp1 = 0;
            tmp4 = 0;
        else
            tmp1 = gvpum*A(vu1,nn);
            tmp4 = gvpum*B(vu1,nn);
        end
        if vu2 < 1 || vm < 1 || um < -vm || um > vm
            tmp2 = 0;
            tmp5 = 0;
        else
            tmp2 = gvmum*A(vu2,nn);
            tmp5 = gvmum*B(vu2,nn);
        end
        if vu3 < 1 || um < -v || um > v
            tmp3 = 0;
            tmp6 = 0;
        else
            tmp3 = gvum*B(vu3,nn);
            tmp6 = gvum*A(vu3,nn);
        end

        A(vu,npnp) = (1/gppnn)*(tmp1 + tmp2 + tmp3);
        B(vu,npnp) = (1/gppnn)*(tmp4 + tmp5 + tmp6);

    end
    end
end

% n+1,-(n+1)
for n=1:(L2-1),
    nn = n^2;
    np = n+1;
    npnm = np^2;
    gpmnn = gpm(n,-n);

    for v = 1:(2*L-n-1),
    for u = (-v):v,      
        vp = v+1;
        vm = v-1;
        up = u+1;    
        vu = v^2 + v + u;
        vu1 = vp^2 + vp + up;
        vu2 = vm^2 + vm + up;
        vu3 = v^2 + v + up;
        gvpup = gmm(vp,up);
        gvmup = gpm(vm,up);
        gvup = gom(v,up);

        if vu1 > tot3 || up < -vp || up > vp
            tmp1 = 0;
            tmp4 = 0;
        else
            tmp1 = gvpup*A(vu1,nn);
            tmp4 = gvpup*B(vu1,nn);
        end
        if vu2 < 1 || vm < 1 || up < -vm || up > vm
            tmp2 = 0;
            tmp5 = 0;
        else
            tmp2 = gvmup*A(vu2,nn);
            tmp5 = gvmup*B(vu2,nn);
        end
        if vu3 > tot3 || up < -v || up > v
            tmp3 = 0;
            tmp6 = 0;
        else
            tmp3 = gvup*B(vu3,nn);
            tmp6 = gvup*A(vu3,nn);
        end

        A(vu,npnm) = (1/gpmnn)*(tmp1 + tmp2 + tmp3);
        B(vu,npnm) = (1/gpmnn)*(tmp4 + tmp5 + tmp6);            

    end
    end
end

% middle
for n=1:(L2-1),
for m=-n:n,
    np = n+1;
    nm = n-1;    
    nn = n^2 + n + m;
    nm = nm^2 + nm + m;
    np = np^2 + np + m;    
    lpnm = lamp(n,m);
    lmnm = lamm(n,m);
    lonm = lamo(n,m);    
    
    for v = 1:(2*L-n-1),
    for u = (-v):v,      
        vp = v+1;
        vm = v-1;    
        vu = v^2 + v + u;
        vu1 = vm^2 + vm + u;
        vu2 = vp^2 + vp + u;    
        lpvmu = lamp(vm,u);
        lmvpu = lamm(vp,u);
        lovu = lamo(v,u);

        if nm < 1 || m < -nm || m > nm
            tmp1 = 0;
            tmp2 = 0;
        else
            tmp1 = -lmnm*A(vu,nm);
            tmp2 = -lmnm*B(vu,nm);
        end
        if vu1 < 1 || vm < 1 || u < -vm || u > vm
            tmp3 = 0;
            tmp4 = 0;
        else
            tmp3 = lpvmu*A(vu1,nn);
            tmp4 = lpvmu*B(vu1,nn);
        end
        if vu2 > tot3 || u < -vp || u > vp
            tmp5 = 0;
            tmp6 = 0;
        else
            tmp5 = lmvpu*A(vu2,nn);
            tmp6 = lmvpu*B(vu2,nn);
        end

        A(vu,np) = (1/lpnm)*(tmp1-lonm*B(vu,nn) ...
                    + tmp3 + tmp5 + lovu*B(vu,nn));
        B(vu,np) = (1/lpnm)*(tmp2-lonm*A(vu,nn) ...
                    + tmp4 + tmp6 + lovu*A(vu,nn));    

    end
    end        
end
end


% resize A and B clipping the extra required for recursion
A = A(1:tot1,:);
B = B(1:tot1,:);

end

function [eta] = epp(n,m)
eta = -bpp(n,m)/(n+1);
end

function [eta] = epm(n,m)
eta = -bpm(n,m)/(n+1);
end

function [eta] = eop(n,m)
eta = gop(n,m);
end

function [eta] = eom(n,m)
eta = gom(n,m);
end

function [eta] = emp(n,m)
eta = bmp(n,m)/n;
end

function [eta] = emm(n,m)
eta = bmm(n,m)/n;
end

function [g] = gpp(n,m)
g = (n/(n+1))*bpp(n,m);
end

function [g] = gpm(n,m)
g = (n/(n+1))*bpm(n,m);
end

function [g] = gop(n,m)
i = complex(0,1);
g = (i/(n*(n+1)))*sqrt((n-m)*(n+m+1));
end

function [g] = gom(n,m)
i = complex(0,1);
g = (i/(n*(n+1)))*sqrt((n+m)*(n-m+1));
end

function [g] = gmp(n,m)
g = ((n+1)/n)*bmp(n,m);
end

function [g] = gmm(n,m)
g = ((n+1)/n)*bmm(n,m);
end

function [z] = zp(n,m)
z = -ap(n,m)/(n+1);
end

function [z] = zo(n,m)
i = complex(0,1);
z = (i*m)/(n*(n+1));
end

function [z] = zm(n,m)
z = am(n,m)/n;
end

function [l] = lamp(n,m)
l = (n/(n+1))*ap(n,m);
end

function [l] = lamo(n,m)
i = complex(0,1);
l = (i*m)/(n*(n+1));
end

function [l] = lamm(n,m)
l = ((n+1)/n)*am(n,m);
end

function [a] = am(n,m)
a = sqrt((n+m)*(n-m)/(2*n+1)/(2*n-1));
end

function [a] = ap(n,m)
a = -sqrt((n+m+1)*(n-m+1)/(2*n+1)/(2*n+3));
end

function [b] = bmm(n,m)
b = -sqrt((n+m)*(n+m-1)/(2*n+1)/(2*n-1));
end

function [b] = bpm(n,m)
b = -sqrt((n-m+2)*(n-m+1)/(2*n+1)/(2*n+3));
end

function [b] = bpp(n,m)
b = sqrt((n+m+2)*(n+m+1)/(2*n+1)/(2*n+3));
end

function [b] = bmp(n,m)
b = sqrt((n-m)*(n-m-1)/(2*n+1)/(2*n-1));
end
