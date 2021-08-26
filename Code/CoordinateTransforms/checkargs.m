function [yn] = checkargs(b,a)
% Coordinate transforms helper function

yn = ~((a == 3 && b == 3) || (a == 6 && b == 3));
if yn 
    disp('not enough input/ouput arguments')
end

