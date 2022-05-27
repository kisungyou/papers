function [center, variation] = corr_mean(input, maxiter, thr)

% CORR_MEAN computes Fréchet mean and variation given an object processed by
% 'corr_initialize' function. Variation is analogous to sample (isotropic)
% variance of data on CORR manifold. 
%
%   * USAGE
%       [center, variation] = corr_mean(input)
%       [center, variation] = corr_mean(input, maxiter)
%       [center, variation] = corr_mean(input, maxiter, thr)
%
%   * INPUT
%       input     an object from 'corr_initialize' for (p,p,N) data.
%       maxiter   (optional) the number of maximum iterations; default 100.
%       thr       (optional) stopping criterion for stopping iterations. 
%                            Default is 1e-8.
%
%   * OUTPUT
%       center    an empirical mean matrix of size (p,p).
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
    error('* corr_mean : incorrect input. Please use CORR_INITIALIZE before usage.');
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
mean_old = mean(input.data, 3);

%% STEP 2. COMPUTE - MEAN
tvecs     = zeros(p,p,N);

for it=1:maxiter
    %   1. compute logarithmic map
    for i=1:N
        tvecs(:,:,i) = corraux_log(mean_old, input.data(:,:,i));
    end
    %   2-2. mean
    mean_tvecs = mean(tvecs,3);
    %   2-3. update using exponential map
    mean_new   = corraux_exp(mean_old, mean_tvecs);
    %   2-4. compute increment
    if (norm(mean_old-mean_new, 'fro') < 1e-10)
        increment = thr/2;
    else
        increment = corraux_dist(mean_old, mean_new);
    end
    %   2-5. update old as new
    mean_old   = mean_new;
    %   2-6. iteration control
    if (increment < thr)
        break;
    end
end
center = mean_old;

%% STEP 3. COMPUTE - FRECHET VARIATION
variation_vec = zeros(1,N);
for i=1:N
    variation_vec(i) = corraux_dist(center, input.data(:,:,i))^2;
end
variation = mean(variation_vec);

end