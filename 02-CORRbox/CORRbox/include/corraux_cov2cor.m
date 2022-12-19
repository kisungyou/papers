function output = corraux_cov2cor(covmat)

A = (covmat + covmat')/2;
invsqrt = diag(1./sqrt(diag(A)));
output = invsqrt*A*invsqrt;

end