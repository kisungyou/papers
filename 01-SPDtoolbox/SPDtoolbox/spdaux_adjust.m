function finout = spdaux_adjust(tgt, thr)
p   = size(tgt,1);
if (nargin<2)
    thr = (p^2)*eps();
end
if (~issymmetric(tgt))
    tgt = (tgt + tgt')/2;
end
[V,D] = eig(tgt);
vecd  = diag(D);
vecd(vecd<thr) = thr;
output = V*diag(vecd)*V';

rk = rank(output);
if (rk>=p)
    finout = output;
else
    finout = spdaux_adjust_recursive(V,D,thr,p);
end
end


function output = spdaux_adjust_recursive(V,D,thr,p)
rk = 0;
while (rk < p)
    thr = p*thr;
    vecd = diag(D);
    vecd(vecd<thr) = thr;
    output = V*diag(vecd)*V';
    rk = (rank(output));
end
end
