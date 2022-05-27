# CORRbox 

**CORRbox** is a `MATLAB` toolbox for learning with correlation matrices (CORR), 
using the Riemannian geometry of correlation manifold. This is performed by viewing 
the space of correlation matrices as a quotient of SPD manifold under the affine-invariant Riemannian metric (AIRM).

### 1. Requirements
This package requires the following MATLAB modules,

  - Optimization toolbox
  - Statistics and Machine Learning Toolbox
  
for fully using the capabilities delivered in the package.

### 2. Add the package
Don't forget to run 
```
addpath(genpath(pwd))
``` 
on MATLAB console when you are navigating this package here.

### 3. Workflow 
When you have CORR-valued data, first it needs to be wrapped as a specific struct that is repeatedly used for all exposed functions in our package, using **corr_initialize** function. Without initialization would functions listed below not work.

### 4. Available functions in `src` folder (in an alphabetical order)
  
| Algorithm | Description |
| :------- | :----------- |
| **`corr_CH`**| a measure of cluster validity by Calinski and Harabasz (1974).|
| **`corr_cmds`**| classical multidimensional scaling. |
| **`corr_eqtestelem`**| element-wise testing for equality of two distributions based on resampling. |
| **`corr_initialize`**| wrap 3d array into a struct after data checking. |
| **`corr_kmeans`**| clustering with $k$-means algorithm. |
| **`corr_kmedoids`**| clustering with $k$-medoids algorithm. |
|**`corr_mean`**| compute Fréchet mean and variation. |
| **`corr_median`**| compute Fréchet median and variation. |
|**`corr_mmds`**| metric multidimensional scaling. |
|**`corr_mse2signal`**| computes mean squared error (MSE) for two given signals.|
|**`corr_pdist`**| compute pairwise distance between every pair of observations.|
|**`corr_pga`**| Principal Geodesic Analysis.|
|**`corr_silhouette`**| a Silhouette measure of cluster validity.|
|**`corr_specc`**| clustering with spectral clustering algorithm.|
|**`corr_test2bg`**| two-sample hypothesis test using Biswas and Ghosh (2014).|
|**`corr_test2wass`**| two-sample hypothesis test using Wasserstein distance.|

  
### 5. Example codes

| Script | Description |
| :------- | :----------- |
| **`example_central`** | compare two measures of central tendency |
| **`example_clustvis`** | clustering and data visualization |
| **`example_eqdist`** | test equality of distributions |
