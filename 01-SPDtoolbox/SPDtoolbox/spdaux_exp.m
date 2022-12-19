function Y = spdaux_exp(X, eta, t)

% exponential mapping

if (nargin < 3)
    t = 1.0;
end
symm = @(x) .5*(x+x');
Y    = symm(X*real(expm(X\(t*eta))));

end