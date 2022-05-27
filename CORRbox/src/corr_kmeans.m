function [idx, C, cost] = corr_kmeans(input, K, initlabel)

% CORR_KMEANS is an adaptation of classical k-means algorithm to
% manifold-valued data on SPD. Lloyd's algorithm was applied in this
% function and it starts using either random initialization or k-means++
% algorithm randomly unless 'initlabel' is provided.
%
%   * USAGE
%       [idx, C, cost] = corr_kmeans(input, K)
%       [idx, C, cost] = corr_kmeans(input, K, initlabel)
%   * INPUT
%       input       an object from 'corr_initialize' for (p,p,N) data.
%       K           predefined number of clusters.
%       initlabel   (optional) length-N vector of class labels.
%   * OUTPUT
%       idx         a length-N vector of class label.
%       C           (p,p,K) 3d array of cluster centers.
%       cost        cost function (weighted sum of per-class variation).
%   * AUTHOR     Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [08/2021] initial implementation.

%% PREPROCESSING
%   1. input
if (~corraux_checker(input))
    error('* corr_kmeans : incorrect input. Please use CORR_INITIALIZE before usage.');
end
p = input.size(1); % dimension
N = input.size(3); % number of data points
%   2. K
K = round(K);
if ((K<=1)||(K>=N))
    error('* corr_kmeans : K should be a number in (1,#{data}).');
end
%   3. initlabel
if (nargin < 3)
    if (rand()<0.5) % 50% chance of random selection
        flagger = true;
        while (flagger)
            initlabel = randi(K,[N,1]);
            flagger = (~(length(unique(initlabel))==K));
        end        
    else            % 50% chance of kmeans++
        flagger = true;
        while (flagger)
            initlabel = corr_kmeans_initialization(input, K);
            flagger   = (~(length(unique(initlabel))==K));
        end
    end
    initlabel = corr_kmeans_initialization(input, K);
else
    initlabel = round(initlabel);
    if (length(initlabel)~=N)
        error('* spd_kmeans : provided INITLABEL should be length of #{data}.');
    end
end


%% Main Computation
%   set up stopping criterion
par_miter = 1000;
par_incre = 0.01;  

%   set 'old' things
old_label = initlabel;                             % length-K vector
[old_means, old_variation] = corr_kmeans_perclass(input, old_label); % oh
old_cost  = corr_kmeans_cost(old_variation, old_label);

% run
itercount = 1;
increment = 10000;
while (increment > par_incre)
    %   1. assignment step give old_means
    pdist2    = corraux_pdist2(old_means, input.data);
    new_label = zeros(size(old_label));
    for i=1:N
        tgtvec = pdist2(:,i);
        tgtidx = find(tgtvec==min(tgtvec));
        if length(tgtidx)>1 % if multiple, choose random one
            new_label(i) = tgtidx(randi(length(tgtidx),1));
        else
            new_label(i) = tgtidx;
        end
    end
    
    %   2. update
    [new_means, new_variation] = corr_kmeans_perclass(input, new_label);
    new_cost = corr_kmeans_cost(new_variation, new_label);
    % Updating is wrong to increase cost function
    % Since it's possible for some runs, only stop it after 25 iterations.
    if ((itercount>=25)&&(new_cost > old_cost)) 
        break;
    end
    increment     = old_cost - new_cost;
    old_means     = new_means;
    old_variation = new_variation;
    old_label     = new_label;
    old_cost      = new_cost;
    
    %   3. stop with itercount
    itercount = itercount + 1;
    if (itercount > par_miter)
        break;
    end
end

%% 3. wrap up for reporting the results
idx  = old_label;
C    = old_means;
cost = old_cost;


end

%% auxiliary functions 
%   (1) initialization
function init_label = corr_kmeans_initialization(input, K)
    p = input.size(1);
    N = input.size(3);
    
    % data reformulation 
    mydata = zeros(N, p*p);
    for i=1:N
        tgtlg       = logm(input.data(:,:,i)); % equivariant embedding
        mydata(i,:) = reshape(tgtlg, [1, p*p]);
    end
    
    % apply matlab's original kmeans function
    init_label = kmeans(mydata,K);
end
%   (2) compute means for each class using 'corraux_mean3d'
function [centers, variation] = corr_kmeans_perclass(input, label)
    ulabel = sort(unique(label));
    K      = length(ulabel);
    p      = size(input.data, 1);
    
    centers   = zeros(p,p,K);
    variation = zeros(K,1);
    
    for i=1:K % main iteration
        idx = (label==ulabel(i));
        [part_mat, part_var] = corraux_mean3d(input.data(:,:,idx));
        centers(:,:,i) = part_mat;
        variation(i)   = part_var;
    end
end
%   (3) compute cost function given variation vector and class label
function cost = corr_kmeans_cost(variation, label)
    ulabel = sort(unique(label));
    K      = length(ulabel);
    
    cost   = 0;
    for k=1:K
        cost = cost + variation(k)*sum(label==ulabel(k));
    end
end