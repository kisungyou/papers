function [label, centers] = corr_kmedoids(input, k)

% CORR_KMEDOIDS performs a k-medoids clustering on the given sample of 
% N observations on the correlation manifold.
%
%   * USAGE
%       [centers, label] = corr_kmedoids(input, k)
%
%   * INPUT
%       input     an object from 'corr_initialize' for (p,p,N) data.
%       k         the number of desired partitions.
%
%   * OUTPUT
%       label     a length-N vector of cluster membership.
%       centers   a (p,p,k) array of cluster centroids stacked as slices.
%
%   * AUTHOR   Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.

%% PREPROCESSING
if (~corraux_checker(input))
    error('* corr_kmedoids : incorrect input. Please use CORR_INITIALIZE before usage.');
end
k = max(1, round(k));
p = input.size(1);
N = input.size(3);

%% COMPUTATION
%  compute the mean
[frechet_mean,~] = corraux_mean3d(input.data);

%  compute the local coordinates
new_data = zeros(N,p*p);
for n=1:N
    mat_tgt = corraux_log(frechet_mean, input.data(:,:,n));
    new_data(n,:) = mat_tgt(:);
end

% apply k-medoids
[idx,~,~,~,midx] = kmedoids(new_data, k);

%% RETURN THE DATA
label   = idx;
centers = input.data(:,:,midx);

end