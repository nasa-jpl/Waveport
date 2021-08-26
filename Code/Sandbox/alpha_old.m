function [alpha] = alpha(L,Lp,k,rji,rgstr)
% 3D scalar translation matrix, alpha

% L,Lp:     maximum order l in row and col
% k:		background wave number
% rji:		[x,y,z] vector from origin to receiving frame
% rgstr:    'rg' for beta = Rg{alpha}, uses spherical Bessel j_l(kr)
%
% alpha:    Translation matrix of size N1xN2
%           N1 = L^2 + 2*L + 1
%           N2 = Lp^2 + 2*Lp + 1
%
% Dependencies: cart2sph, sbesselh, sbesselj

% convert rji to spherical coordinates
[r th phi] = cart2sph(rji(1),rji(2),rji(3));

% Determine maximum order and transpose condition
flip = 0;
if Lp > L
    tmp = Lp;
    Lp = L;
    L = tmp;
    flip = 1;
end
L3 = L + Lp;
tot1 = L^2 + 2*L + 1;
tot2 = Lp^2 + 2*Lp + 1;
tot3 = L3^2 + 2*L3 + 1;
alpha = zeros(tot3,tot2);
str = 'mono';

% alpha_{vu,00}
ylm = sphericalY(L3,th,phi,'mono');
if nargin == 5 && strcmp(rgstr,'rg')
    bess = sbesselj(0:L3,k*r);
else
    bess = sbesselh(0:L3,k*r);
end
for v=0:L3,
for u=-v:v,
    vu = v^2 + v + u + 1;
    alpha(vu,1) = sqrt(4*pi)*(-1)^(v)*conj(ylm(vu))*bess(v+1);
end
end

% b_nn+ alpha_vu,n+1,n+1 = ...
for n = 0:(Lp-1),
    nn = n^2 + 2*n + 1;
    np = n+1;        
    npnp = np^2 + 2*np + 1;
    bpnn = bp(n,n);    
    for v = 0:(L3-1-n),
    for u = (-v):v,
        in1 = v^2 + v + u + 1;
        vp = v+1;
        vm = v-1;
        um = u-1;
        vu1 = vm^2 + vm + um + 1;
        vu2 = vp^2 + vp + um + 1;
        bp1 = bp(vm,um);
        bm1 = bm(vp,um);        
        if vm < 0 || um < -vm
            alpha(in1,npnp) = (1/bpnn)*(bm1*alpha(vu2,nn));
        else
            alpha(in1,npnp) = (1/bpnn)*(bp1*alpha(vu1,nn) + bm1*alpha(vu2,nn));
        end    
    end
    end
end

% b_n,-n+ alpha_vu,n+1,-(n+1) = ...
for n = 0:(Lp-1),
    nn = n^2 + 1;
    np = n+1;
    npnp = np^2 + 1;
    bpnn = bpm(n,-n);    
    for v = 0:(L3-1-n),
    for u = (-v):v,
        in1 = v^2 + v + u + 1;       
        vp = v+1;
        vm = v-1;
        up = u+1;
        vu1 = vm^2 + vm + up + 1;
        vu2 = vp^2 + vp + up + 1;
        bp1 = bpm(vm,up);
        bm1 = bmm(vp,up);
        if vm < 0 
            alpha(in1,npnp) = (1/bpnn)*(bm1*alpha(vu2,nn));
        else
            alpha(in1,npnp) = (1/bpnn)*(bp1*alpha(vu1,nn) + bm1*alpha(vu2,nn));
        end    
    end
    end
end


% a_nm+ alpha_vu,n+1,m = -a_nm- alpha_vu,n1,m ...
for n = 0:(Lp-1),
for m = -n:n,
    np = n+1;
    nm = n-1;
    nm0 = n^2 + n + m + 1;
    nmm = nm^2 + nm + m + 1;
    npm = np^2 + np + m + 1;    
    ap1 = ap(n,m);
    am1 = am(n,m);        
    for v = 0:(L3-n-1),
    for u = (-v):v,          
        vp = v+1;
        vm = v-1;
        vu0 = v^2 + v + u + 1;
        vu1 = vm^2 + vm + u + 1;
        vu2 = vp^2 + vp + u + 1;
        ap2 = ap(vm,u);
        am2 = am(vp,u);
        if (m >= 0 && nm < 0) || (m<0 && nm <= 0)
            if vm < 0 || u < -vm
                alpha(vu0,npm) = (1/ap1)*(am2*alpha(vu2,nm0));
            else
                alpha(vu0,npm) = (1/ap1)*(ap2*alpha(vu1,nm0) + am2*alpha(vu2,nm0));
            end    
        else        
            if vm < 0 || u < -vm
                alpha(vu0,npm) = (1/ap1)*(-am1*alpha(vu0,nmm)+am2*alpha(vu2,nm0));
            else
                alpha(vu0,npm) = (1/ap1)*(-am1*alpha(vu0,nmm)+ap2*alpha(vu1,nm0) + am2*alpha(vu2,nm0));
            end
        end
    end
    end
end
end

% 
% % a_nm+ alpha_vu,n+1,m = -a_nm- alpha_vu,n1,m ...
% for n = 0:(Lp-1),
% for m = 0:n,
%     np = n+1;
%     nm = n-1;
%     nm0 = n^2 + n + m + 1;
%     nmm = nm^2 + nm + m + 1;
%     npm = np^2 + np + m + 1;    
%     ap1 = ap(n,m);
%     am1 = am(n,m);        
%     for v = 0:(L3-n-1),
%     for u = (-v):v,          
%         vp = v+1;
%         vm = v-1;
%         vu0 = v^2 + v + u + 1;
%         vu1 = vm^2 + vm + u + 1;
%         vu2 = vp^2 + vp + u + 1;
%         ap2 = ap(vm,u);
%         am2 = am(vp,u);
%         if nm < 0    
%             if vm < 0 || u < -vm
%                 alpha(vu0,npm) = (1/ap1)*(am2*alpha(vu2,nm0));
%             else
%                 alpha(vu0,npm) = (1/ap1)*(ap2*alpha(vu1,nm0) + am2*alpha(vu2,nm0));
%             end    
%         else        
%             if vm < 0 || u < -vm
%                 alpha(vu0,npm) = (1/ap1)*(-am1*alpha(vu0,nmm)+am2*alpha(vu2,nm0));
%             else
%                 alpha(vu0,npm) = (1/ap1)*(-am1*alpha(vu0,nmm)+ap2*alpha(vu1,nm0) + am2*alpha(vu2,nm0));
%             end
%         end
%     end
%     end
% end
% end
% 
% % a_nm+ alpha_vu,n+1,m = -a_nm- alpha_vu,n1,m ...
% for n = 1:(Lp-1),           
% for m = (-n):(-1),    
%     np = n+1;
%     nm = n-1;    
%     nm0 = n^2 + n + m + 1;
%     nmm = nm^2 + nm + m + 1;
%     npm = np^2 + np + m + 1;    
%     ap1 = ap(n,m);
%     am1 = am(n,m);        
%     for v = 0:(L3-n-1),
%     for u = (-v):v,          
%         vp = v+1;
%         vm = v-1;        
%         vu0 = v^2 + v + u + 1;
%         vu1 = vm^2 + vm + u + 1;
%         vu2 = vp^2 + vp + u + 1;    
%         ap2 = ap(vm,u);
%         am2 = am(vp,u);
%         if nm <= 0    
%             if vm < 0 || u < -vm
%                 alpha(vu0,npm) = (1/ap1)*(am2*alpha(vu2,nm0));
%             else
%                 alpha(vu0,npm) = (1/ap1)*(ap2*alpha(vu1,nm0) + am2*alpha(vu2,nm0));
%             end    
%         else        
%             if vm < 0 || u < -vm
%                 alpha(vu0,npm) = (1/ap1)*(-am1*alpha(vu0,nmm)+am2*alpha(vu2,nm0));
%             else
%                 alpha(vu0,npm) = (1/ap1)*(-am1*alpha(vu0,nmm)+ap2*alpha(vu1,nm0) + am2*alpha(vu2,nm0));
%             end    
%         end    
%     end
%     end
% end
% end

% trim the matrix
alpha = alpha(1:tot1,:);

% apply (-1)^(n+v+m+u) if flip
if flip
    alpha2 = zeros(tot1,tot2);
    for n=0:Lp,
    for m=-n:n,
        nm = n^2 + n + m + 1;
        nmm = n^2 + n - m + 1;
        for v=0:L,
        for u=-v:v,
            uv = v^2 + v + u + 1;
            uvm = v^2 + v - u + 1;
            alpha2(uv,nm) = (-1)^(n+v+m+u)*alpha(uvm,nmm);
        end
        end
    end
    end
	alpha = alpha2.';
end

end


function amc = am(n,m)
    amc = sqrt((n+m)*(n-m)/(2*n+1)/(2*n-1));
end

function apc = ap(n,m)
    apc = -sqrt((n+m+1)*(n-m+1)/(2*n+1)/(2*n+3));
end

function bmc = bm(n,m)
    bmc = sqrt((n-m)*(n-m-1)/(2*n+1)/(2*n-1));
end

function bmc = bmm(n,m)
    bmc = -sqrt((n+m)*(n+m-1)/(2*n+1)/(2*n-1));
end

function bpc = bp(n,m)
    bpc = sqrt((n+m+2)*(n+m+1)/(2*n+1)/(2*n+3));
end

function bpc = bpm(n,m)
    bpc = -sqrt((n-m+2)*(n-m+1)/(2*n+1)/(2*n+3));
end