function ind = box2ind(l,i)
% Column index for 1D fast multipole method expansion coefficients
%
% l,i:  level and subdivision number
%
% ind:  column index

ind = 2^l - 2 + i ;


