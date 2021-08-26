
% 
% L = 10;
% tot= L^2 + 2*L;
% 
% Q = zeros(2*tot,2*tot);
% rgQ = zeros(2*tot,2*tot);
% 
% indM = 1:tot;
% indN = (tot+1):(2*tot);
% 
% [tab] = lmtable(L);
% l = tab(:,1);
% l2 = 1:L;
% 
% 
% for n=1:N,
% 
% % convert surface point to spherical coordinates
% [rp thp php] = cart2sph(X(n),Y(n),Z(n));  
% 
% % unconjugated vector spherical harmonics
% [V2, U2, ~, ~] = BC(L,thp,php,'norm');
% W2 = sphericalY(L,thp,php);
%     
% [c1 p1 b1] = ebcmbesselfunc(tot,L,l,l2,k1,r,[]);
% [rgc1 rgp1 rgb1] = ebcmbesselfunc(tot,L,l,l2,k1,r,'rg');
% 
% 
% end

% 
% % loop over surface points (slow but saves memory)
% for n=1:N,
%     if mod(n,1000) == 0
%         mydisp(n,N)
%     end
%     % convert surface point to spherical coordinates
%     [rp thp php] = cart2sph(X(n),Y(n),Z(n));  
%     
%     % convert spherical unit vectors at this point to Cartesian components
%     [rhat_x rhat_y rhat_z]      = sph2cart(rp,thp,php,1,0,0);
% 	[thhat_x thhat_y thhat_z]   = sph2cart(rp,thp,php,0,1,0);
%     [phhat_x phhat_y phhat_z]   = sph2cart(rp,thp,php,0,0,1);
% 
%     % dot product of surface normal and spherical unit vectors
%     % include the differential surface element here
%     n_dot_r  = dS(n)*dot(n_hat(n,:),[rhat_x rhat_y rhat_z],2);
% 	n_dot_th = dS(n)*dot(n_hat(n,:),[thhat_x thhat_y thhat_z],2);
%     n_dot_ph = dS(n)*dot(n_hat(n,:),[phhat_x phhat_y phhat_z],2);
% 
%     % unconjugated vector spherical harmonics (row arrays)
%     [V2, U2, ~, ~] = BC(L,thp,php,'norm');
%     W2 = sphericalY(L,thp,php);
% 
%     % conjugated vector spherical harmonics (row arrays)
%     V1 = conj(V2);
%     U1 = conj(U2);
%     W1 = conj(W2);
%     
%     % bessel functions (row arrays)
%     [c1 p1 b1]          = ebcmbesselfunc(tot,L,l,l2,k1,rp,[]);
%     [rgc1 rgp1 rgb1]    = ebcmbesselfunc(tot,L,l,l2,k1,rp,'rg');
%     [rgc2 rgp2 rgb2]    = ebcmbesselfunc(tot,L,l,l2,k2,rp,'rg');
%     
%     % cross product pre-multiplications (row arrays)
%     [cU1 cV1 bU1 bV1 pW1]           = ebcmprod(c1,p1,b1,U1,V1,W1);
%     [rgcU1 rgcV1 rgbU1 rgbV1 rgpW1] = ebcmprod(rgc1,rgp1,rgb1,U1,V1,W1);
%     [rgcU2 rgcV2 rgbU2 rgbV2 rgpW2] = ebcmprod(rgc2,rgp2,rgb2,U2,V2,W2);
%     
%     %[cU1 cV1 bU1 bV1 pW1]           = ebcmprod(1,1,1,U1,V1,W1);
%     %[rgcU1 rgcV1 rgbU1 rgbV1 rgpW1] = ebcmprod(1,1,1,U1,V1,W1);
%     %[rgcU2 rgcV2 rgbU2 rgbV2 rgpW2] = ebcmprod(1,1,1,U2,V2,W2);
%         
%     % Outer products for Q: (l,m) index along rows, (p,q) index along columns
%     M1hat_cross_rgM2 =  n_dot_r*(cV1.'*rgcU2 - cU1.'*rgcV2);
%     N1hat_cross_rgN2 =  n_dot_r*(bV1.'*rgbU2 - bU1.'*rgbV2) + n_dot_th*(bU1.'*rgpW2 - pW1.'*rgbU2) + n_dot_ph*(pW1.'*rgbV2 - bV1.'*rgpW2);
%     M1hat_cross_rgN2 =  n_dot_r*(cU1.'*rgbU2 + cV1.'*rgbV2) - n_dot_th*(cV1.'*rgpW2) - n_dot_ph*(cU1.'*rgpW2);
%     N1hat_cross_rgM2 = -n_dot_r*(bU1.'*rgcU2 + bV1.'*rgcV2) + n_dot_th*(pW1.'*rgcV2) + n_dot_ph*(pW1.'*rgcU2);
%     
%     % Q
% 	Q(indM,indM) = Q(indM,indM) + const*M1hat_cross_rgN2 + N1hat_cross_rgM2; % need + to match sphere T-matrix
% 	Q(indM,indN) = Q(indM,indN) + const*M1hat_cross_rgM2 + N1hat_cross_rgN2;
%     Q(indN,indM) = Q(indN,indM) + const*N1hat_cross_rgN2 + M1hat_cross_rgM2;
%     Q(indN,indN) = Q(indN,indN) + const*N1hat_cross_rgM2 + M1hat_cross_rgN2;
%     
%     % Outer products for RgQ
% 	rgM1hat_cross_rgM2 =  n_dot_r*(rgcV1.'*rgcU2 - rgcU1.'*rgcV2);
% 	rgN1hat_cross_rgN2 =  n_dot_r*(rgbV1.'*rgbU2 - rgbU1.'*rgbV2) + n_dot_th*(rgbU1.'*rgpW2 - rgpW1.'*rgbU2) + n_dot_ph*(rgpW1.'*rgbV2 - rgbV1.'*rgpW2);
% 	rgM1hat_cross_rgN2 =  n_dot_r*(rgcU1.'*rgbU2 + rgcV1.'*rgbV2) - n_dot_th*(rgcV1.'*rgpW2) - n_dot_ph*(rgcU1.'*rgpW2);
% 	rgN1hat_cross_rgM2 = -n_dot_r*(rgbU1.'*rgcU2 + rgbV1.'*rgcV2) + n_dot_th*(rgpW1.'*rgcV2) + n_dot_ph*(rgpW1.'*rgcU2);
%    
%     % RgQ
%     RgQ(indM,indM) = RgQ(indM,indM) + const*rgM1hat_cross_rgN2 + rgN1hat_cross_rgM2;
%     RgQ(indM,indN) = RgQ(indM,indN) + const*rgM1hat_cross_rgM2 + rgN1hat_cross_rgN2;
%     RgQ(indN,indM) = RgQ(indN,indM) + const*rgN1hat_cross_rgN2 + rgM1hat_cross_rgM2;
%     RgQ(indN,indN) = RgQ(indN,indN) + const*rgN1hat_cross_rgM2 + rgM1hat_cross_rgN2;
%     
%  
% end



% ebcm
[Tmm Tmn Tnm Tnn] = ebcm(L,X,Y,Z,dS,nx,ny,nz,k1,k2);


