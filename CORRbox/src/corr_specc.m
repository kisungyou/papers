function label = corr_specc(input, k, nnbd)

% CORR_SPECC applies a spectral clustering algorithm by a version of 
% Zelnik-Manor and Perona (2005) where an affinity matrix is constructed 
% using the local distance to the 'nnbd'-th nearest observations for N 
% observations on the correlation manifold.
%
%   * USAGE
%       label = corr_specc(input, k, nnbd)
%
%   * INPUT
%       input     an object from 'corr_initialize' for (p,p,N) data.
%       k         the number of desired partitions.
%       nnbd      (Optional) the neighborhood size.
%
%   * OUTPUT
%       label     a length-N vector of cluster membership.
%
%   * AUTHOR   Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.


%% PREPROCESSING
if (~corraux_checker(input))
    error('* corr_specc : incorrect input. Please use CORR_INITIALIZE before usage.');
end
if (nargin < 3)
    nnbd = 7;
end
k = max(1, round(k));
N = input.size(3);
nnbd = min(max(round(nnbd), 1), N-1);

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

%   compute the nearest distance vector
vec_nearest = zeros(1,N);
for n=1:N
    tgt_sorted     = sort(mat_pdist(n,:));
    vec_nearest(n) = tgt_sorted(nnbd);
end

%   compute the similarity matrix using Zelnik-Manor
mat_S = ones(N,N);
for i=1:(N-1)
    for j=(i+1):N
        mat_S(i,j) = exp(-(mat_pdist(i,j)^2)/(vec_nearest(i)*vec_nearest(j)));
        mat_S(j,i) = mat_S(i,j);
    end
end

%   perform spectral clustering using NJW
label = spectralcluster(mat_S, k, 'Distance','precomputed','LaplacianNormalization','symmetric');

end