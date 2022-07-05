function H = spdaux_log(X, Y)

% logarithmic map
symm = @(x) .5*(x+x');
H    = symm(X*real(logm(X\Y)));

end
