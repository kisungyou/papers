1. Notes
    - Add this folder by 'addpath(genpath(pwd))'
    - Every function (not in src folder) has documentation. 
    - All SPD matrices should be stacked as 3d array of size (p,p,N).
    - For the data, use "spd_initialize" to wrap it securely.

2. Available Functions
    - spd_eqtest     : testing equality of two distributions based on permutation testing.
    - spd_initialize : wrap 3d array into a struct after data checking.
    - spd_kmeans     : perform k-means clustering on SPD manifold.
    - spd_mean       : compute Fréchet mean and variation.
    - spd_pdist      : compute pairwise distance between data pairs.
    - spd_pga        : Principal Geodesic Analysis.
    - spd_smooth     : apply kernel smoothing on SPD matrices.
    - spd_smoothcv   : cross-validation for bandwidth selection of 'spd_smooth'.
    - spd_validity   : numerical scores to determine which 'clustering' is good.
    
