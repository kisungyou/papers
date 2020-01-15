function score = spd_clustval(input, label, name)

% SPD_CLUSTVAL provides a measure of validating cluster quality defined
% accordingly with respect to each measure. Default validity measure is
% 'Silhouette' index by Rousseeux (1987). 
%
%   * USAGE
%       score = spd_clustval(input, label)
%       score = spd_clustval(input, label, name)
%   * INPUT
%       input       an object from 'spd_initialize' for (p,p,N) data.
%       label       a length-N vector of class labels, or
%                   a cell containing multiple length-N label vectors.
%       name        (optional) name of validity index. Following measures
%                   are currently supported. (+) means the higher a score
%                   is, the better label is.
%           'Silhouette' (+) Silhouette index
%   * OUTPUT
%       score       a constant ('label' is vector), or
%                   a vector   ('label' is cell of vectors) of scores.
%   * AUTHOR     Kisung You (kyou@nd.edu)
%   * HISTORY
%       0.1. [01/2019] initial implementation.
%		0.2. [08/2019] change the function name for clarification.



%% Preprocessing : Check Inputs
%   0. Default as Silhouette
if (nargin < 3)
    name = 'Silhouette';
end
%   1. input
if (~spdaux_checker(input))
    error('* spd_clustval : incorrect input. Please use SPD_INITIALIZE before usage.');
end
N = input.size(3);
%   2. label
if (isvector(label)&&(~iscell(label)))
    labelclass = 'single';
    if (length(label)~=N)
        error('* spd_clustval : LABEL should be a vector of length : input.size(3).');
    end
elseif (isvector(label)&&iscell(label))
    labelclass = 'multiple';
    ncell = length(label);
    for i=1:ncell
        if (~isvector(label{i}))
            error('* spd_clustval : each element in a cell LABEL should be a vector.');
        end
        if (length(label{i})~=N)
            error('* spd_clustval : each element in a cell LABEL should be a vector of length : input.size(3).');
        end
    end
else
    error('* spd_clustval : LABEL should be either a vector or a cell of vectors.');
end

%% Main Computation
if strcmp(labelclass,'single')
    score = spd_validity_single(input, round(label), name);
else
    score = zeros(length(label),1);
    for i=1:length(label)
        score(i) = spd_validity_single(input, round(label{i}), name);
    end
end
end


%% Branching Function
function score = spd_validity_single(input, lvector, name)
    myname = lower(name);
    switch myname
        case 'silhouette'
            score = spd_validity_MAXSil(input, lvector);
        otherwise
            error(strcat('* spd_clustval : '',name,'' is not supported yet.'));
    end
end

%% (max) 1. Silhouette
function score = spd_validity_MAXSil(input, label)
    ulabel = sort(unique(label));
    K      = length(ulabel);
    if (K<2)
        error('* spd_clustval : ''Silhouette'' index does not work for K=1.');
    end
    n      = input.size(3);
    
    vec_a   = zeros(n,1); % a(i)
    vec_b   = zeros(n,1); % b(i)
    
    indexer = spdaux_indexer(label); % index per cluster
    distmat = spd_pdist(input);      % pairwise distance matrix
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
