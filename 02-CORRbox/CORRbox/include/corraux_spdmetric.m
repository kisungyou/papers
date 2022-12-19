function mval = corraux_spdmetric(X, eta1, eta2)
    sol1 = (X\eta1);
    sol2 = (X\eta2);
    mval = sum(diag(sol1'*sol2)); 
end