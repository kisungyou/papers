function distmat = spd_pdist(input,mode)

% SPD_PDIST returns a symmetric matrix where each (i,j)-th element is a
% distance measure between i- and j-th elements from given SPD data.
%
%   * USAGE
%       distmat = spd_pdist(input, mode)
%   * INPUT
%       input       an object from 'spd_initialize' for (p,p,N) data.
%       mode        either 'intrinsic' or 'extrinsic'.
%   * OUTPUT
%       distmat     an (N,N) matrix recording distance data pairs.
%   * AUTHOR     Kisung You (kyou@nd.edu)
%   * HISTORY
%       0.1. [01/2019] initial implementation.

%% Preprocessing : Check Inputs
if (~spdaux_checker(input))
    error('* spd_pdist : incorrect input. Please use SPD_INITIALIZE before usage.');
end
N = input.size(3);

if (nargin < 2)
    mymode = 'intrinsic';
else
    mymode = mode;
end

%% Main Computation
distmat = zeros(N,N);
for i=1:(N-1)
    tgt1 = input.data(:,:,i);
    for j=(i+1):N
        tgt2 = input.data(:,:,j);
        dval = spdaux_dist(tgt1,tgt2,mymode);
        
        distmat(i,j) = dval;
        distmat(j,i) = dval;
    end
end
end
