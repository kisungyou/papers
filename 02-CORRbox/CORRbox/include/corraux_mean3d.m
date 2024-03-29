function [center, variation] = corraux_mean3d(array3d)

% simply compute mean matrix given 3d array

%% preprocessing
thr = 1e-8;
p   = size(array3d,1);
N   = size(array3d,3);

%% step 1 : extrinsic mean in the ambient space

mean_old = mean(array3d, 3);

%% step 2 : compute intrinsic mean from extrinsic one

increment = 10000;
tvecs     = zeros(p,p,N);
iter      = 1;
while (increment > thr)
    %   2-1. compute logarithmic map
    for i=1:N
        tvecs(:,:,i) = corraux_log(mean_old, array3d(:,:,i));
    end
    %   2-2. mean
    mean_tvecs = mean(tvecs,3);
    %   2-3. update using exponential map
    mean_new   = corraux_exp(mean_old, mean_tvecs);
    %   2-4. compute increment
    increment  = corraux_dist(mean_old, mean_new);
    %   2-5. update old as new
    mean_old   = mean_new;
    %   2-6. iteration count
    iter = iter + 1;
    if (iter >= 1000)
        break;
    end
end

%% step 3. compute Frechet variation
center        = mean_old;
variation_vec = zeros(1,N);
for i=1:N
    variation_vec(i) = corraux_dist(center, array3d(:,:,i))^2;
end
variation     = mean(variation_vec);

end