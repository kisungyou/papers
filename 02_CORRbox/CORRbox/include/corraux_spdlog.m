function H = corraux_spdlog(X,Y)
    symm = @(x) .5*(x+x');
    H    = symm(X*real(logm(X\Y)));
end