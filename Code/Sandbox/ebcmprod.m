function [cU cV bU bV pW] = ebcmprod(c,p,b,U,V,W);
    
cU = c.*U;
cV = c.*V;
bU = b.*U;
bV = b.*V;
pW = p.*W;