function outval = corraux_metric(C, X)

term1 = corraux_spdmetric(C,X,X);

n = size(C,1);
Cinv = inv(C);
AC = C.*Cinv;

mu = ((eye(n)+AC)\diag(diag(C\X)))*ones(n,1);

term2 = -2*mu'*(eye(n)+AC)*mu;
outval = term1+term2;

end

% function mval = corraux_metric(X, eta1, eta2)
%     sol1 = (X\eta1);
%     sol2 = (X\eta2);
%     mval = sum(diag(sol1'*sol2)); 
% end