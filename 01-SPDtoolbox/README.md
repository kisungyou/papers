# SPDtoolbox

### 1. Download

In order to download this folder only, please use the following commands.
```
git clone --depth 1 --filter=blob:none --no-checkout https://github.com/kisungyou/papers
cd papers/
git checkout master -- 01-SPDtoolbox
```

### 2. Add the package
Don't forget to run 
```
addpath(genpath(pwd))
``` 
on MATLAB console when you are navigating this package here.


### 3. Workflow 
When you have SPD-valued data, first it needs to be wrapped as a specific struct that is repeatedly used for all exposed functions in our package, using **spd_initialize** function. Without initialization would functions listed below not work.

### 4. Choice of Metric : AIRM and LERM
Our package supports two types of metrics

  * AIRM : Affine-Invariant Riemannian Metric, and
  * LERM : Log-Euclidean Riemannian Metric
  
in the sense that AIRM is an **_intrinsic_** geometry while LERM is governed by **_extrinsic_** geometry on the manifold of symmetric and positive-definite matrices.

### 5. Available functions in `src` folder (in an alphabetical order)

| Algorithm | Description |
| :------- | :----------- |
|**`spd_clustval`**| a measure of validating cluster quality|
|**`spd_eqtest`**| testing equality of two distributions based on permutation testing.|
|**`spd_eqtestelem`**| element-wise testing equality of two distributions based on resampling.|
|**`spd_icca`**| extracts independent covariance components.|
|**`spd_initialize`**| wrap 3d array into a struct after data checking.|
|**`spd_kmeans`**| perform $k$-means clustering on SPD manifold.|
|**`spd_mean`**| compute Fréchet mean and variation.|
|**`spd_mse2signal`**| computes mean squared error (MSE) for two given signals.|
|**`spd_pdist`**| compute pairwise distance between data pairs.|
|**`spd_pga`**| Principal Geodesic Analysis.|
|**`spd_smooth`**| apply kernel smoothing on SPD matrices.|
|**`spd_smoothcv`**| cross-validation for bandwidth selection of 'spd_smooth'.|
  
### 6. Example codes
| Script | Description |
| :------- | :----------- |
|**`example_eqdist.m`** | test equality of distributions under both AIRM and LERM.|
|**`example_mean`**   | compare intrinsic and extrinsic Fréchet means.|
|**`example_kmeans`** | applying k-means algorithm and cluster validity index.|
