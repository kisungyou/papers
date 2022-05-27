function [embedding, stress] = corr_mmds(input, ndim)

% CORR_MMDS performs a metric multidimensional scaling for visualizing 
% the distribution of N observation on the correlation manifold using the 
% Kruskal's normalized stress criterion.
%
%   * USAGE
%       [embedding,stress] = corr_mmds(input, ndim)
%
%   * INPUT
%       input     an object from 'corr_initialize' for (p,p,N) data.
%       ndim      the target dimension.
%
%   * OUTPUT
%       embedding an (N,ndim) embedding coordinates.
%       stress    the normalized stress value at the end of iterations.
%
%   * AUTHOR   Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.

%% PREPROCESSING
if (~corraux_checker(input))
    error('* corr_mmds : incorrect input. Please use CORR_INITIALIZE before usage.');
end
if (nargin < 2)
    ndim = 2;
end
ndim = max(1, round(ndim));
N    = input.size(3);

%% COMPUTATION
%   compute the pairwise distance matrix
mat_pdist = zeros(N,N);
for i=1:(N-1)
    tgt1 = input.data(:,:,i);
    for j=(i+1):N
        tgt2 = input.data(:,:,j);
        dval = corraux_dist(tgt1,tgt2);
        
        mat_pdist(i,j) = dval;
        mat_pdist(j,i) = dval;
    end
end

% apply 'cmdscale' function from Statistics & Machine Learning Toolbox
[embedding, stress] = mdscale(mat_pdist, ndim);

end
