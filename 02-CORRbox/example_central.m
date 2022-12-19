%% EXAMPLE : Fréchet mean and median
%  data preparation
%   (p,p,N) array where slices are correlation matrix of sample data 
%   from normal distribution. Change the first part for different settings.

%% Setting and data generation
clear; close all; clc; addpath(genpath(pwd));

p = 10;  % parameter : dimension
N = 10;  % parameter : number of slices

array3d = zeros(p,p,N);
for i=1:10
    array3d(:,:,i) = corr(randn(10*p,p)*2);
end
corr3d = corr_initialize(array3d); % Don't forget to initialize

%% Compute measures of central tendency
[fmean, var1]   = corr_mean(corr3d);    % Fréchet mean
[fmedian, var2] = corr_median(corr3d);  % Fréchet median


%% Visualize
%  title text
title_var1 = sprintf("Fréchet mean with variation=%.4f",var1);
title_var2 = sprintf("Fréchet median with variation=%.4f",var2);

%  draw
figure()
subplot(1,2,1); imagesc(fmean); colorbar; title(title_var1); axis square;
subplot(1,2,2); imagesc(fmedian); colorbar; title(title_var2); axis square;

