function distmat = corr_pdist(input)

% CORR_PDIST returns a symmetric matrix where each (i,j)-th element is a
% distance measure between i- and j-th elements from given CORR data.
%
%   * USAGE
%       distmat = corr_pdist(input)
%
%   * INPUT
%       input       an object from 'corr_initialize' for (p,p,N) data.
%
%   * OUTPUT
%       distmat     an (N,N) matrix recording distance data pairs.
%
%   * AUTHOR     Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.

%% Preprocessing : Check Inputs
if (~corraux_checker(input))
    error('* corr_pdist : incorrect input. Please use CORR_INITIALIZE before usage.');
end
N = input.size(3);

%% Main Computation
distmat = zeros(N,N);
for i=1:(N-1)
    tgt1 = input.data(:,:,i);
    for j=(i+1):N
        tgt2 = input.data(:,:,j);
        dval = corraux_dist(tgt1,tgt2);
        
        distmat(i,j) = dval;
        distmat(j,i) = dval;
    end
end

end