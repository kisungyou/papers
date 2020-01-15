%% EXAMPLE : k-means clustering
%% data preparation
%   Generate 20 covariance matrices.
%   First 10 slices are covariance of sample data from normal distribution.
%   Latter 10 slices are those from uniform distribution.

array3d = zeros(5,5,20);
for i=1:10
    array3d(:,:,i) = corr(randn(200,5)*2);
end
for j=11:20
    array3d(:,:,j) = corr(rand(200,5)) + diag(ones(5,1));
end

spd3d = spd_initialize(array3d); % Don't forget to initialize

%% Example 1 : Compute Pairwise Distance
distmat = spd_pdist(spd3d);
imagesc(distmat); colorbar; % visualize

%% Example 2 : compute mean and variation
%  Since data are from 2 distinct ones, variation should be large.
[center, variation] = spd_mean(spd3d);

%% Example 3 : k-means
[idx, C, cost] = spd_kmeans(spd3d, 2);

%% Example 4 : validity
%   4-1. compute a score given a single label
score = spd_clustval(spd3d, idx);

%   4-2. compare multiple clustering labels' quality
idx2to6 = cell(5,1);
for i=1:5
    [tmpidx,~,~] = spd_kmeans(spd3d, i+1); % iterate K=2 to K=6
    idx2to6{i} = tmpidx;
end
score5 = spd_clustval(spd3d, idx2to6);

%% Visualize
figure(1)
subplot(1,2,1); imagesc(distmat); colorbar; title("1. Pairwise Distance Matrix"); axis square;
subplot(1,2,2); plot(2:6, score5, "-o"); title("2. Silhouette Index (k=2 should be highest)"); xlabel("k : number of clusters"); axis square;
set(gcf, "color","white");