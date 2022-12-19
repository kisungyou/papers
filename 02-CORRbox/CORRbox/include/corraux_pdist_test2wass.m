function pvalue = corraux_pdist_test2wass(DXY, m, n, niter)

if (nargin < 4)
    nreps = 999;
else
    nreps = max(9, round(niter));
end

%% COMPUTE THE STATISTIC

m = round(m);
n = round(n);


%% COMPUTATION - STATISTIC

DZ0 = DXY(1:m, (m+1):(m+n));
Tmn = corraux_wassersteinD(DZ0, 2.0); 

%% COMPUTATION - ITERATIONS
printer = true;
Tvec = zeros(1,nreps);
for i=1:nreps
    % random permutation
    id0 = randperm(m+n); % random permutation
    idx = id0(1:m);
    idy = id0((m+1):(m+n));
    
    DZ1 = DXY(idx,idy);
    Tvec(i) = corraux_wassersteinD(DZ1, 2.0);
    if (printer)
        fprintf("* test2wass : iteration %d/%d complete.. \n",i,nreps);
    end
end

%% P-VALUE
pvalue = (sum(Tvec>Tmn)+1)/(nreps+1);
end