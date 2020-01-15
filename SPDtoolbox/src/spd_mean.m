function [center, variation] = spd_mean(input, thr, mode)

% SPD_MEAN computes Fréchet mean and variation given an object processed by
% 'spd_initialize' function. Variation is analogous to sample (isotropic)
% variance of data on SPD manifold. To be more specific, this function
% computes 'intrinsic' mean and corresponding variation.
%
%   * USAGE
%       [center, variation] = spd_mean(input)
%       [center, variation] = spd_mean(input, thr)
%       [center, variation] = spd_mean(input, thr, mode)
%   * INPUT
%       input     an object from 'spd_initialize' for (p,p,N) data.
%       thr       (optional) stopping criterion for stopping iterations. 
%                            Default is 1e-8.
%       mode      type of analysis. 'extrinsic' or 'intrinsic'.
%   * OUTPUT
%       center    an empirical mean matrix of size (p,p).
%       variation sample Fréchet variation.
%
%   * AUTHOR   Kisung You (kyou@nd.edu)
%   * HISTORY
%       0.1. [01/2019] initial implementation.
%
%   See also SPD_INITIALIZE

%% checker
if (~spdaux_checker(input))
    error('* spd_mean : incorrect input. Please use SPD_INITIALIZE before usage.');
end

%% thr
if (nargin < 2)
    thr = 1e-8;
end
if (nargin < 3)
    mode = 'extrinsic';
end

%% step 1. compute extrinsic mean as an initializer
p = input.size(1); % dimension
N = input.size(3); % number of slices

array3d_logm = zeros(p,p,N);
for i=1:N
    array3d_logm(:,:,i) = real(logm(input.data(:,:,i)));
end
mean_old = real(expm(mean(array3d_logm,3)));
mean_old = spdaux_adjust(mean_old);

%% step 2. compute intrinsic mean from extrinsic one

if (strcmp(mode,'extrinsic'))
    center = mean_old;
elseif (strcmp(mode,'intrinsic'))
    increment = 10000;
    tvecs     = zeros(p,p,N);
    iter      = 1;
    while (increment > thr)
        %   2-1. compute logarithmic map
        for i=1:N
            tvecs(:,:,i) = spdaux_log(mean_old, input.data(:,:,i));
        end
        %   2-2. mean
        mean_tvecs = mean(tvecs,3);
        %   2-3. update using exponential map
        mean_new   = spdaux_exp(mean_old, mean_tvecs);
        %   2-4. compute increment
        if (norm(mean_old-mean_new, 'fro') < 1e-6)
            increment = thr/2;
        else
            increment = spdaux_dist(mean_old, mean_new, mode);
        end
        %   2-5. update old as new
        mean_old   = mean_new;
        %   2-6. iteration count
        iter = iter + 1;
        if (iter >= 1000)
            break;
        end
    end
    center = mean_old;
else
    error('* spd_mean : MODE should be either INTRINSIC or EXTRINSIC.');
end
%% step 3. compute Frechet variation
variation_vec = zeros(1,N);
for i=1:N
    variation_vec(i) = spdaux_dist(center, input.data(:,:,i), mode)^2;
end
variation     = mean(variation_vec);

end
