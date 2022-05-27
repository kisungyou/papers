function [center, variation] = corr_median(input, maxiter, thr)

% CORR_MEDIAN computes Fréchet median and variation given an object processed by
% 'corr_initialize' function. Variation is analogous to sample (isotropic)
% variance of data on CORR manifold. 
%
%   * USAGE
%       [center, variation] = corr_median(input)
%       [center, variation] = corr_median(input, maxiter)
%       [center, variation] = corr_median(input, maxiter, thr)
%
%   * INPUT
%       input     an object from 'corr_initialize' for (p,p,N) data.
%       maxiter   (optional) the number of maximum iterations; default 100.
%       thr       (optional) stopping criterion for stopping iterations. 
%                            Default is 1e-8.
%
%   * OUTPUT
%       center    an empirical median matrix of size (p,p).
%       variation sample Fréchet variation.
%
%   * AUTHOR   Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.
%
%   See also CORR_INITIALIZE


%% PREPROCESS
%  check the input
if (~corraux_checker(input))
    error('* corr_median : incorrect input. Please use CORR_INITIALIZE before usage.');
end
%  threshold value
if (nargin < 2)
    maxiter = 100;
end
if (nargin < 4)
    thr = 1e-8;
end
maxiter = max(50, round(maxiter));

%% STEP 1. INITIALIZE
p = input.size(1);
N = input.size(3);
median_old = corraux_pseudoext(input.data);

%% STEP 2. COMPUTE - MEDIAN
vec_log   = zeros(p,p,N);
vec_dist  = zeros(1,N);

for it=1:maxiter
    % 1. compute log and distance
    for n=1:N
        [tmpLogmat, tmpdist] = corr_median_dual(median_old, input.data(:,:,n));
        vec_log(:,:,n) = tmpLogmat;
        vec_dist(n)    = tmpdist;
    end
    % 2. updating information
    tmpLog = zeros(p,p);
    tmpNum = 0;
    for n=1:N
        if (vec_dist(n) > 1e-10)
            tmpNum = tmpNum + (1/vec_dist(n));
            tmpLog = tmpLog + (vec_log(:,:,n)/vec_dist(n));
        end
    end
    % 3. do the update
    median_new = corraux_exp(median_old, tmpLog/tmpNum, 1.0);
    increment  = norm(median_old-median_new,"fro");
    median_old = median_new;
    % 4. stopping criterion & update
    if (increment < thr)
        break;
    end
end
center = median_old;

%% STEP 3. COMPUTE - FRECHET VARIATION
variation_vec = zeros(1,N);
for i=1:N
    variation_vec(i) = corraux_dist(center, input.data(:,:,i));
end
variation = mean(variation_vec);

end


%% compute log and distance at the same time
function [logXY, dval] = corr_median_dual(X,Y)
    D = corraux_findD(X, Y);
    Z = D*Y*D;
    logXY = corraux_log(X, Z);
    dval  = corraux_spddist(X, Z);
end