0. Add the package
   Don't forget to run
   >> addpath(genpath(pwd))
   when you are navigating this package here.

1. Workflow 
   When you have SPD-valued data, first it needs to be wrapped as a specific struct that is repeatedly used for all exposed functions in our package.

2. Choice of Metric : AIRM and LERM
   Our package supports two types of metrics
	- AIRM : Affine-Invariant Riemannian Metric, and
	- LERM : Log-Euclidean Riemannian Metric
   in the sense that AIRM is an 'intrinsic' geometry while LERM is governed by 'extrinsic' geometry on the manifold of symmetric and positive-definite matrices.

3. Available functions in 'src' folder (in an alphabetical order)
    - spd_clustval   : a measure of validating cluster quality
    - spd_eqtest     : testing equality of two distributions based on permutation testing.
    - spd_eqtestelem : elementwise testing equality of two distributions based on resampling.
    - spd_icca       : extracts independent covariance components.
    - spd_initialize : wrap 3d array into a struct after data checking.
    - spd_kmeans     : perform k-means clustering on SPD manifold.
    - spd_mean       : compute Fréchet mean and variation.
    - spd_mse2signal : computes mean squared error (MSE) for two given signals.
    - spd_pdist      : compute pairwise distance between data pairs.
    - spd_pga        : Principal Geodesic Analysis.
    - spd_smooth     : apply kernel smoothing on SPD matrices.
    - spd_smoothcv   : cross-validation for bandwidth selection of 'spd_smooth'.

4. Example codes
    - example_eqdist : test equality of distributions under both AIRM and LERM.
    - example_mean   : compare intrinsic and extrinsic Frechet means
    - example_kmeans : applying k-means algorithm and cluster validity index.
    
