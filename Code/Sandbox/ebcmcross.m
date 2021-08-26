function [M1xN2 M1xM2 N1xN2 N1xM2] = ebcmcross(cU1,cV1,bU1,bV1,pW1,cU2,cV2,bU2,bV2,pW2,n_dot_r,n_dot_th,n_dot_ph);
    
M1xN2 =  n_dot_r*(cU1.'*bU2 + cV1.'*bV2) - n_dot_th*(cV1.'*pW2) - n_dot_ph*(cU1.'*pW2);
M1xM2 =  n_dot_r*(cV1.'*cU2 - cU1.'*cV2);
N1xN2 =  n_dot_r*(bV1.'*bU2 - bU1.'*bV2) + n_dot_th*(bU1.'*pW2 - pW1.'*bU2) + n_dot_ph*(pW1.'*bV2 - bV1.'*pW2);
N1xM2 = -n_dot_r*(bU1.'*cU2 + bV1.'*cV2) + n_dot_th*(pW1.'*cV2) + n_dot_ph*(pW1.'*cU2);
