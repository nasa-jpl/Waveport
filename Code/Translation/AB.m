function [A B] = AB(L,Lp,k,rji,rgstr,normstr)
% Vector translation matrices via the scalar translation matrix
%
% L,Lp:     Maximum row/column harmonic degree 
% k:        Background wave number
% rji:      [x,y,z] vector from the source to receiver frame
% rgstr:    'rg' for the regular form: Rg{A} and Rg{B}
% normstr:  'norm' for fully normalized wave functions
%
% A,B:      Translation matrices size N1 x N2
%           N1 = L^2 + 2*L
%           N2 = Lp^2 + 2*Lp
%
% Dependencies: lm2ind, alpha

if nargin == 4
    rgstr = [];
    normstr = [];
end
% convert rji to spherical coordinates
[r th phi] = cart2sph(rji(1),rji(2),rji(3));

% dimensions of the matrix
tot1 = L^2 + 2*L;  % row dimension
tot2 = Lp^2 + 2*Lp;  % col dimension

% preallocation
A = zeros(tot1,tot2);
B = zeros(tot1,tot2);

% compute alpha (or beta = Rg{alpha}) for L1+1, and L2+1
[alphaLp1] = alpha(L+1,Lp+1,k,rji,rgstr);

% the columns of A and B are shifted one column to the right 
% relative to the columns of alpha.
grab = 2:(tot2+1);

% constants used in the interations
const1 = k*r*sin(th)*exp(-i*phi)/2;
const2 = -const1;
const3 = -k*r*sin(th)*exp(i*phi)/2;
const4 = -const3;
const5 = k*r*cos(th);
const6 = const5;
const7 = 1i*const5;
const8 = 1i*k*r*sin(th)*exp(i*phi)/2;
const9 = 1i*k*r*sin(th)*exp(-i*phi)/2;

% loop over rows of the vector matrices
for v=1:L,
for u=-v:v,   
    % precomputed linear indecies
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
    vu8 = v^2 + v + up + 1;
    vu9 = v^2 + v + um + 1;
    
    % iteration coefficients
    c1 = const1*(1/(v+1))*sqrt((v-u+2)*(v-u+1)/(2*v+1)/(2*v+3));
    c2 = const2*(1/v)*sqrt((v+u-1)*(v+u)/(2*v-1)/(2*v+1));
    c3 = const3*(1/(v+1))*sqrt((v+u+2)*(v+u+1)/(2*v+1)/(2*v+3));
    c4 = const4*(1/v)*sqrt((v-u)*(v-u-1)/(2*v-1)/(2*v+1));    
    c5 = const5*(1/(v+1))*sqrt((v+u+1)*(v-u+1)/(2*v+1)/(2*v+3));
    c6 = const6*(1/v)*sqrt((v+u)*(v-u)/(2*v-1)/(2*v+1));    
    c7 = const7*u*(1/(v*(v+1)));
    c8 = const8*(1/(v*(v+1)))*sqrt((v-u)*(v+u+1));
    c9 = const9*(1/(v*(v+1)))*sqrt((v+u)*(v-u+1));
    
    % conditional statements to handle index limits
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
    if up < -v || up > v, tmp8 = 0;
    else tmp8 = c8*alphaLp1(vu8,grab); end
    if um < -v || um > v, tmp9 = 0; 
    else tmp9 = c9*alphaLp1(vu9,grab); end
            
    % translation matricies
    A(vu-1,:) = alphaLp1(vu,grab) + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6;
    B(vu-1,:) = c7*alphaLp1(vu,grab) + tmp8 + tmp9;
end
end

% apply normalization if needed
if nargin == 6 && strcmp(normstr,'norm')
    for n=1:Lp,
    for mp=-n:n,
    for l=1:L,
    for m=-l:l,
        const = sqrt(l*(l+1))/sqrt(n*(n+1));
        A(lm2ind(l,m),lm2ind(n,mp)) = const*A(lm2ind(l,m),lm2ind(n,mp));
        B(lm2ind(l,m),lm2ind(n,mp)) = const*B(lm2ind(l,m),lm2ind(n,mp));
    end
    end
    end
    end
end


