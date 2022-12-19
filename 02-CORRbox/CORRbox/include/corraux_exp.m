function expmat = corraux_exp(C, X, t)

% ray starts at C in the direction of X

if (nargin < 3)
    t = 1.0;
end

horX     = corraux_proj_horizontal(C, X);
Chalf    = sqrtm(C);
Chalfinv = inv(Chalf);
 
expmat = corraux_cov2cor(Chalf*expm(t*Chalfinv*horX*Chalfinv)*Chalf);


end


% function Y = corraux_exp(X, eta, t)
% 
% % exponential mapping
% 
% if (nargin < 3)
%     t = 1.0;
% end
% symm = @(x) .5*(x+x');
% Ytmp = symm(X*real(expm(X\(t*eta))));
% Y    = corraux_cov2cor(Ytmp);
% 
% end