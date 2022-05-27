function score = corr_CH(input, label)

% corr_CH is a measure of clustering quality that was proposed by 
% Calinski and Harabasz in 1974. We interpret the score that the higher the score is, 
% the better the given clustering is.
%
%   * USAGE
%       score = corr_CH(input, label)
%
%   * INPUT
%       input       an object from 'corr_initialize' for (p,p,N) data.
%       label       a length-N vector of class labels.
%
%   * OUTPUT
%       score       a CH score.
%
%   * AUTHOR     Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [09/2021] initial implementation.

%% PREPROCESSING
%   1. input
if (~corraux_checker(input))
    error('* corr_CH : incorrect input. Please use CORR_INITIALIZE before usage.');
end
p = input.size(1);
N = input.size(3);
%   2. label
if (~isvector(label))
    error("* corr_CH : input LABEL should be a vector.");
end
if (length(label)~=N)
    error('* corr_CH : LABEL should be a vector of length : input.size(3).');
end
label = round(label);

%% MAIN COMPUTATION
%  Preliminary
ulabel = sort(unique(label));
K      = length(ulabel);
if (K<2)
    error('* corr_CH : index does not work for K=1.');
end

%  Compute Means
mean_group = zeros(p,p,K);
for k=1:K
    idk = find(label==ulabel(k));
    [mean_tmp,~] = corraux_mean3d(input.data(:,:,idk));
    mean_group(:,:,k) = mean_tmp;
end
[mean_all,~]   = corraux_mean3d(input.data);

% Compute Separation (numerator)
val_numerator = 0;
for k=1:K
    val_numerator = val_numerator + (corraux_dist(mean_all, mean_group(:,:,k))^2)*sum(label==ulabel(k))/(K-1);
end

% Compute Cohesion   (denominator)
val_denominator = 0;
for k=1:K
    tmp_sum  = 0;
    tmp_id   = find(label==ulabel(k));
    tmp_nobs = length(tmp_id);
    for it=1:tmp_nobs
        tmp_sum = tmp_sum + (corraux_dist(mean_group(:,:,k), input.data(:,:,tmp_id(it)))^2)/(N-K);
    end
    val_denominator = val_denominator + tmp_sum;
end

% return the score 
score = val_numerator/val_denominator;

end

% https://www.geeksforgeeks.org/calinski-harabasz-index-cluster-validity-indices-set-3/
