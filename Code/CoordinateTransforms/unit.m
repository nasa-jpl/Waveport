function varargout = unit(varargin)
% Unit vector
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
% varargout:    unit vector components matching the inputs

if length(varargin) == 1 
    vec = varargin{:};
    Ndim = length(vec(1,:));
    m = mag(vec);
    vec = vec./repmat(m,1,Ndim);
    vec((m==0),:) = 0; 
    varargout{1} = vec;
elseif length(varargin)>1
    m = mag(varargin{:});
    Ndim = length(varargin)
    for n=1:Ndim,
       varargout{n} =  varargin{n}./m;
    end
end

