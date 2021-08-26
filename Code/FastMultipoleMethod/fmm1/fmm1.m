function [f] = fmm1(ak,S)
% 1D fast multipole method on x = [-1 1] for kernel of the form  
%
%       f(x_j) = sum_k^N a_k/(x_j - x_k)
%
% Specialized for non-overlapping source/observation points.  Takes 
% inputs from the preparatory funciton fmm1prep.
%
% ak:       source amplitudes (real or complex)
% S:        precomputed structure from fmm1prep
%
% f:        function evaluted at observation points xj
%
% dependencies: fmm1prep, fmm1u, box2ind

f = zeros(S.Nj,1);
S.Phi(:) = 0;
S.Psi(:) = 0;
ak = ak(:);

% compute far-field at nlevs boxes, k points
const = 2^S.nlevs - 2;
for b = 1:length(S.xk_notempty),
    i = S.xk_notempty(b);
    ind = const + i;   	% ind = 2^S.nlevs - 2 + i;  box2ind(S.nlevs,i)
    S.Phi(:,ind) = S.uk_far{i}*ak(S.xk_ind{i});
end
 
% far-field at each subinterval, each level
for l=(S.nlevs-1):-1:2,
    tl = 2^l;
    tlp1m3 = 2^(l+1) - 3;
   for i=1:tl,
       ind1 = tl - 2 + i;    % ind1 = 2^l - 2 + i;        box2ind(l,i)
       ind2 = tlp1m3 + 2*i;  % ind2 = 2^(l+1) + 2*i - 3;  box2ind(l+1,2*i-1)
       ind3 = ind2 + 1;      % ind3 = ind2 + 1;           box2ind(l+1,2*i)
       S.Phi(:,ind1) = S.ML*S.Phi(:,ind2) + S.MR*S.Phi(:,ind3);
   end    
end

% local expansions each level, each subinterval
for l=1:(S.nlevs-1),
    tl = 2^l;
    tlp1 = 2^(l+1);
    for i=1:(2^l),
        ti = 2*i;
        tli = tlp1 + ti;
        if (ti - 3) >= 1	% (2*i-3) >= 1
            ind = tli - 5;  % ind = 2^(l+1) + 2*i - 5; box2ind(l+1,2*i-3)
            tmp1L = S.T2*S.Phi(:,ind);
            tmp1R = S.T1*S.Phi(:,ind);
        else
            tmp1L = 0;
            tmp1R = 0;
        end
        if (ti-2) >= 1      % (2*i-2) >= 1
            ind = tli - 4;  % ind = 2^(l+1) + 2*i - 4; box2ind(l+1,2*i-2);
            tmp2 = S.T2*S.Phi(:,ind);
        else
            tmp2 = 0;
        end
        if (ti+1) <= tlp1 	% (2*i+1) <= 2^(l+1)
            ind = tli - 1;  % ind = 2^(l+1) + 2*i - 1; box2ind(l+1,2*i+1)
            tmp3 = S.T3*S.Phi(:,ind);
        else
            tmp3 = 0;
        end
        if (ti+2) <= tlp1	% (2*i+2) <= 2^(l+1)
            ind = tli;      % ind = 2^(l+1) + 2*i; box2ind(l+1,2*i+2)
            tmp4L = S.T4*S.Phi(:,ind);
            tmp4R = S.T3*S.Phi(:,ind);
        else
            tmp4L = 0;
            tmp4R = 0;
        end
        ind1 = tli - 3;     % ind1 = 2^(l+1) + 2*i - 3; box2ind(l+1,2*i-1)
        ind2 = tli - 2;     % ind2 = 2^(l+1) + 2*i - 2; box2ind(l+1,2*i)
        ind3 = tl + i - 2;	% ind3 = 2^l + i - 2;       box2ind(l,i)
        S.Psi(:,ind1)   = S.SL*S.Psi(:,ind3) + tmp1L + tmp3 + tmp4L;
        S.Psi(:,ind2)   = S.SR*S.Psi(:,ind3) + tmp1R + tmp2 + tmp4R;     
     end
end

% evaluate local expansions at subset of xj
const = 2^S.nlevs - 2;
for b = 1:length(S.xj_notempty),
    i = S.xj_notempty(b);
    ind = const + i;	% ind = 2^S.nlevs - 2; box2ind(S.nlevs,i)
    f(S.xj_ind{i}) = S.uj_local{i}*S.Psi(:,ind);
end

% evaluate near-neighbors directly with sparse matrix-vector multiply
f = f + S.M*ak;

end



