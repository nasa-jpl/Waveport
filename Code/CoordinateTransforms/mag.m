function [m] = mag(varargin)
% Magnitude of vector
% 
% varargin:     length(varargin) == 1, 
%               E.g., unit(vec)
%               vec is an [N x Ndim] array of vector components. 
%               N points for Ndim dimensional vector
%
%               length(vararin) > 1
%               E.g., unit(X,Y,Z,W,...)  
%               X,Y,Z,W vector components of any size
%
% m:            vector magintude for each point N

if length(varargin) == 1 
    m = sqrt(sum(varargin{1}.^2,2));
elseif length(varargin)>1
    Ndim = length(varargin);
    m = 0;
    for n=1:Ndim,
       v = varargin{n};
       m = m + v.^2; 
    end
    m = sqrt(m);
end
