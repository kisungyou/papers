function pvalue = corraux_pdist_test2bg(DXY, m, n, niter)

if (nargin < 4)
    nreps = 999;
else
    nreps = max(9, round(niter));
end

%% COMPUTE THE STATISTIC

m = round(m);
n = round(n);

%% COMPUTATION - STATISTIC

id_m = 1:m;
id_n = (m+1):(m+n);

DX0 = DXY(id_m, id_m);
DY0 = DXY(id_n, id_n);
DZ0 = DXY(id_m, id_n);
Tmn = corraux_test2bg_statistic(DX0, DY0, DZ0); 

%% COMPUTATION - ITERATIONS
Tvec = zeros(1,nreps);
for i=1:nreps
    % random permutation
    id0 = randperm(m+n); % random permutation
    
    id_m = id0(1:m);
    id_n = id0((m+1):(m+n));
    
    DX1 = DXY(id_m, id_m);
    DY1 = DXY(id_n, id_n);
    DZ1 = DXY(id_m, id_n);
    Tvec(i) = corraux_test2bg_statistic(DX1,DY1,DZ1);
end

%% P-VALUE
pvalue = (sum(Tvec>=Tmn)+1)/(nreps+1);
end

