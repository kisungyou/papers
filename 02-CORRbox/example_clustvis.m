%% EXAMPLE : clustering and data visualization
%  This script loads the derived model correlation matrices that were used 
%  in our SPD and Correlation paper. We demonstrate visualization as well
%  as clustering results on the perturbed data.

%% Setup and load the model matrices

clear; close all; clc; addpath(genpath(pwd));
load("cluster_three_corrs.mat");

%% Data generation

ncopy = 10;   % number of samples per class
sd    = 0.1;  % degree of perturbation

data_3d = zeros(5,5,3*ncopy);
for i=1:ncopy
    data_3d(:,:,i) = corraux_perturb(C1, sd);
    data_3d(:,:,i+ncopy) = corraux_perturb(C2, sd);
    data_3d(:,:,i+(2*ncopy)) = corraux_perturb(C3, sd);
end
data_obj = corr_initialize(data_3d);
data_lab = repelem(1:3, ncopy);

%% Compute an embedding with metric multidimensional scaling
%  other visualization/dimension reduction methods may be used similarly.

[embed2,~] = corr_mmds(data_obj, 2); disp("* complete 1/4 : metric mds.");

%% Cluster analysis using spectral clustering algorithm

mynbd  = 5;
specc2 = corr_specc(data_obj, 2, mynbd); disp("* complete 2/4 : spectral clustering with K=2.");
specc3 = corr_specc(data_obj, 3, mynbd); disp("* complete 3/4 : spectral clustering with K=3.");
specc4 = corr_specc(data_obj, 4, mynbd); disp("* complete 4/4 : spectral clustering with K=4.");

%% Visualize

tiledlayout(1,3, "TileSpacing", "compact", "Padding", "tight");
nexttile
scatter(embed2(:,1), embed2(:,2), 16, specc2, 'filled');
xlabel("MMDS1"); ylabel("MMDS2"); axis square; title("specc K=2");
nexttile
scatter(embed2(:,1), embed2(:,2), 16, specc3, 'filled');
xlabel("MMDS1"); ylabel("MMDS2"); axis square; title("specc K=3");
nexttile
scatter(embed2(:,1), embed2(:,2), 16, specc4, 'filled');
xlabel("MMDS1"); ylabel("MMDS2"); axis square; title("specc K=4");