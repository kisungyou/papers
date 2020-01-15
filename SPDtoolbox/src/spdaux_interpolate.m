function C3 = spdaux_interpolate(C1,C2,w)

% Description
%   : it computes C3 = w*C1 + (1-w)*C2

symm = @(x) .5*(x+x');
Ct   = symm(C1*real(logm(C1\C2)));
C3   = symm(C1*real(expm(C1\(w*Ct))));


end