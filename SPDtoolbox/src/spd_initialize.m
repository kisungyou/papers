function output = spd_initialize(array3d, thr)

% SPD_INITIALIZE converts a 3d array into a specific struct that is
% repeatedly used for all exposed functions in our package. If any slice of
% provided 3d array is either asymmetric or not positive definite, this
% function automatically 
%
%   * USAGE
%       output = spd_initialize(array3d)
%       output = spd_initialize(array3d, thr)
%   * INPUT
%       array3d     a 3d array of size (p,p,N) where each slice is SPD
%       thr         threshold value to zero out non-positive eigenvalues
%   * OUTPUT
%       output      a struct containing following elements,
%         .name     name of manifold 'spd'
%         .size     a vector of data size
%         .data     (p,p,N) 3d array of processed/checked data
%   * AUTHOR     Kisung You (kyou@nd.edu)
%   * HISTORY
%       0.1. [01/2019] initial implementation.
%       0.2. [01/2019] change to adapt threshold level for eigenvalues.

%% main run
%  1. checkers
%   1-1. 3d array
size3d = size(array3d);
if (length(size3d)~=3)
    error('* spd_initialize : an input needs to be a 3-dimensional array.');
end
%   1-2. no need for single slice
if (size3d(3)==1)
    error('* spd_initialize : single sliced input is not considered.');
end
%   1-3. square matrices
if (size3d(1)~=size3d(2))
    error('* spd_initialize : each slice must be squared matrices');
end
p = size3d(1);
N = size3d(3);



%   3. let's iterate over slices
if (nargin < 2)
    thr = (p^2)*eps();
end
if (thr < eps())
    error('* spd_initialize : THR must be a small, positive number.');
end
mydata = zeros(p,p,N);
for i=1:N
    tgt = array3d(:,:,i);
    if (~issymmetric(tgt))
        tgt = (tgt + tgt')/2;
    end
    tmp = spdaux_adjust(tgt, thr);
    mydata(:,:,i) = (tmp + tmp')/2;
end

%   2. prepare a struct
output = struct();
output.name = 'spd';
output.size = [p,p,N];
output.data = mydata;

end
