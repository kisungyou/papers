function dval = corraux_spddist(X,Y)
trnorm = @(A) sqrt(trace(A*A));
dval   = real(trnorm(real(logm(X\Y)))); % AIRM/LERM; not clear what to use.
end