function [A B] = AB_old(L1,L2,k,rs)
% 3D vector addition theorem (via scalar alpha)

% A and B expand outgoing spherical vector wave function from the origin
%  as incoming waves about rs.  Implementation of Chew's recursive 
% algorithm built on the scalar translation matrix alpha.

% Dependencies: ALPHA.m

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
theta = atan2(sqrt(x.^2+y.^2),z); 
ct = cos(theta);
phi = atan2(y,x);

tot1 = L1^2 + 2*L1;  % row dimension
tot2 = L2^2 + 2*L2;  % col dimension

A = zeros(tot1,tot2);
B = zeros(tot1,tot2);

% need alpha for L1+1, and L2+1
[alphaLp1] = alpha(L1+1,L2+1,k,rs);

% final columns of A and B are shifted in place 
% one column in alpha, this grabs them
grab = 2:(tot2+1);

% constants used below
coef1 = k*r*sin(theta)*exp(-i*phi);
coef2 = -coef1;
coef3 = -conj(coef1);
coef4 = conj(coef1);
coef5 = k*r*cos(theta);
coef81 = i*k*r*sin(theta)*exp(i*phi);
coef82 = i*k*r*sin(theta)*exp(-i*phi);

for v=1:L1,
for u=-v:v,   
    vp = v+1;
    vm = v-1;
    up = u+1;
    um = u-1;          
    vu = v^2 + v + u + 1;
    vu1 = vp^2 + vp + um + 1;
    vu2 = vm^2 + vm + um + 1;
    vu3 = vp^2 + vp + up + 1;
    vu4 = vm^2 + vm + up + 1;
    vu5 = vp^2 + vp + u + 1;
    vu6 = vm^2 + vm + u + 1;
    vu81 = v^2 + v + up + 1;
    vu82 = v^2 + v + um + 1;
    
    c1 = coef1*(1/(2*(v+1)))*...
        sqrt((v-u+2)*(v-u+1)/(2*v+1)/(2*v+3));
    c2 = coef2*(1/(2*v))*...
        sqrt((v+u-1)*(v+u)/(2*v-1)/(2*v+1));
    c3 = coef3*(1/(2*(v+1)))*...
        sqrt((v+u+2)*(v+u+1)/(2*v+1)/(2*v+3));
    c4 = coef4*(1/(2*v))*...
        sqrt((v-u)*(v-u-1)/(2*v-1)/(2*v+1));    
    c5 = coef5*(1/(v+1))*...
        sqrt((v+u+1)*(v-u+1)/(2*v+1)/(2*v+3));
    c6 = coef5*(1/v)*...
        sqrt((v+u)*(v-u)/(2*v-1)/(2*v+1));    
    c7 = coef5*i*u*(1/(v*(v+1)));
    c81 = coef81*(1/(2*v*(v+1)))*sqrt((v-u)*(v+u+1));
    c82 = coef82*(1/(2*v*(v+1)))*sqrt((v+u)*(v-u+1));
    
    
    if um < -vp || um > vp, tmp1 = 0;
    else tmp1 = c1*alphaLp1(vu1,grab); end
    if um < -vm || um > vm, tmp2 = 0;
    else tmp2 = c2*alphaLp1(vu2,grab); end    
    if up < -vp || up > vp, tmp3 = 0;
    else tmp3 = c3*alphaLp1(vu3,grab); end
    if up < -vm || up > vm, tmp4 = 0;
    else tmp4 = c4*alphaLp1(vu4,grab); end    
    if up < -vp || up > vp, tmp5 = 0;
    else tmp5 = c5*alphaLp1(vu5,grab); end       
    if u < -vm || u > vm, tmp6 = 0;
    else tmp6 = c6*alphaLp1(vu6,grab); end   
    if up < -v || up > v, tmp81 = 0;
    else tmp81 = c81*alphaLp1(vu81,grab); end
    if um < -v || um > v, tmp82 = 0; 
    else tmp82 = c82*alphaLp1(vu82,grab); end
            
    A(vu-1,:) = alphaLp1(vu,grab) + tmp1 + tmp2 + tmp3 ...
                + tmp4 + tmp5 + tmp6;
    B(vu-1,:) = c7*alphaLp1(vu,grab) + tmp81 + tmp82;
    
end
end















