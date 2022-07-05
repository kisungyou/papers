function distmat = spdaux_pdist2(array1, array2, mode)

% pairwise distance for (array1) and (array2)

if (nargin < 3)
    mymode = 'extrinsic';
else
    mymode = mode;
end

M = size(array1,3);
N = size(array2,3);

distmat = zeros(M,N);
for i=1:M
    tgt1 = array1(:,:,i);
    for j=1:N
        tgt2 = array2(:,:,j);
        distmat(i,j) = spdaux_dist(tgt1,tgt2,mymode);
    end
end
end
