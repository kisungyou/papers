%% EXAMPLE : Frechet mean
%  data preparation
%   (p,p,N) array where slices are covariance matrix of sample data 
%   from normal distribution. Change the first part for different settings.

clear; close all; clc;
p = 10;   % parameter : dimension
N = 10;  % parameter : number of slices

array3d = zeros(p,p,N);
for i=1:10
    array3d(:,:,i) = cov(randn(10*p,p)*2);
end
spd3d = spd_initialize(array3d); % Don't forget to initialize

%% Compute AIRM / Intrinsic Mean and Variation
[ctd_airm, var_airm] = spd_mean(spd3d, 1e-8, "intrinsic");

%% Compute LERM / Extrinsic Mean and Variation
[ctd_lerm, var_lerm] = spd_mean(spd3d, 1e-8, "extrinsic");

%% Visualize
title_airm = sprintf("AIRM with variation=%.4f",var_airm);
title_lerm = sprintf("LERM with variation=%.4f",var_lerm);

figure()
subplot(1,2,1); imagesc(ctd_airm); colorbar; title(title_airm); axis square;
subplot(1,2,2); imagesc(ctd_lerm); colorbar; title(title_lerm); axis square;
set(gcf, "color","white");