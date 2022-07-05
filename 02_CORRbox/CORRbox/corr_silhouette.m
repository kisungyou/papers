function score = corr_silhouette(input, label)

% CORR_SILHOUETTE provides a measure of validating cluster quality defined
% by 'Silhouette' index by Rousseeux (1987). The higher the score is, the 
% better the given clustering is.
%
%   * USAGE
%       score = corr_silhouette(input, label)
%
%   * INPUT
%       input       an object from 'corr_initialize' for (p,p,N) data.
%       label       a length-N vector of class labels.
%
%   * OUTPUT
%       score       a Silhouette score.
%
%   * AUTHOR     Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.

%% PREPROCESSING
%   1. input
if (~corraux_checker(input))
    error('* corr_silhouette : incorrect input. Please use CORR_INITIALIZE before usage.');
end
N = input.size(3);
%   2. label
if (~isvector(label))
    error("* corr_silhouette : input LABEL should be a vector.");
end
if (length(label)~=N)
    error('* corr_silhouette : LABEL should be a vector of length : input.size(3).');
end
label = round(label);

%% MAIN COMPUTATION
ulabel = sort(unique(label));
K      = length(ulabel);
if (K<2)
    error('* corr_silhouette : index does not work for K=1.');
end

n       = input.size(3);
vec_a   = zeros(n,1); % a(i)
vec_b   = zeros(n,1); % b(i)

indexer = corraux_indexer(label); % index per cluster
distmat = corr_pdist(input);      % pairwise distance matrix
for i=1:n
    % a(i) : average within cluster
    %        singleton index should be taken care.
    idx1     = indexer{(ulabel==label(i))};
    if (length(idx1)==1)
        vec_a(i) = 0;
    else
        vec_a(i) = mean(distmat(i,setdiff(idx1,i)));
    end
    % b(i) : weird.. but gotta do.
    otherlabels = setdiff(ulabel, label(i)); % non-overlapping classes
    mindistance = zeros(K-1,1);
    for k=1:(K-1)
        idx2 = indexer{(ulabel==otherlabels(k))};
        mindistance(k) = mean(distmat(i,idx2));
    end
    vec_b(i) = min(mindistance);
end

% final computation
vec_s = zeros(n,1);
for i=1:n
    term_top = vec_b(i)-vec_a(i);
    term_bot = max(vec_a(i),vec_b(i));
    vec_s(i) = term_top/term_bot;
end
score = mean(vec_s);

end