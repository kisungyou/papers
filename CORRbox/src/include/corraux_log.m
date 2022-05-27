function logmat = corraux_log(C1, C2)

Dtmp  = corraux_findD(C1, C2);
C2tmp = Dtmp*C2*Dtmp;

Sigma = C1;
V     = corraux_spdlog(C1, C2tmp);

Delta_Sigma_inv  = diag(1./sqrt(diag(Sigma)));
Delta_Sigma_inv2 = diag(1./diag(Sigma));


inners = Delta_Sigma_inv2*diag(diag(V))*Sigma + Sigma*diag(diag(V))*Delta_Sigma_inv2;
logmat = Delta_Sigma_inv*(V - 0.5*(inners))*Delta_Sigma_inv;


end



% function logXY = corraux_log(X, Y)
%     D = corraux_findD(X, Y);
%     Z = D*Y*D;
%     logXY = tmp_log(X, Z);
% end

%% temporary function
function H = tmp_log(X,Y)
    symm = @(x) .5*(x+x');
    H    = symm(X*real(logm(X\Y)));
end