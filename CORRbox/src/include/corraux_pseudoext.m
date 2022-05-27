function output = corraux_pseudoext(array3d)

% size argument
p = size(array3d, 1);
N = size(array3d, 3);

% equivariant mapping
array3d_logm = zeros(p,p,N);
for i=1:N
    array3d_logm(:,:,i) = real(logm(array3d(:,:,i)));
end

% compute the output
output = corraux_cov2cor(real(expm(mean(array3d_logm,3))));

end