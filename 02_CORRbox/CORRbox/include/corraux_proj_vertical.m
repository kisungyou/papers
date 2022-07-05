function verV = corraux_proj_vertical(S, V)

% S in SPD(n)
% V in T_S SPD(n)

n    = size(V,1);
Sinv = inv(S);

mu   = (eye(n) + (S.*Sinv))\diag(diag(Sinv*V))*ones(n,1);
phi  = (mu*ones(1,n) + ones(n,1)*mu');
verV = S.*phi;

end