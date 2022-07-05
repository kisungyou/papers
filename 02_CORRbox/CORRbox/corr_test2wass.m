function pvalue = corr_test2wass(input1, input2, niter)

% CORR_TEST2WASS performs a hypothesis testing for the equality of two
% distributions given samples 'input1' and input2' according to the 
% nonparametric method using 2-Wasserstein distance.
%
%   * USAGE
%       pvalue = corr_test2wass(input1, input2)
%       pvalue = corr_test2wass(input1, input2, niter)
%
%   * INPUT
%       input1     an object from 'spd_initialize' for class 1
%       input2     an object from 'spd_initialize' for class 2
%       niter      (optional) the number of operations; default=999.
%
%   * OUTPUT
%       pvalue    p-value under 'H0 : two are equally distributed'.
%
%   * AUTHOR   Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.

%% PREPROCESSING
%  number of permutations
if (nargin < 3)
    nreps = 999;
else
    nreps = max(9, round(niter));
end

%  check inputs
if (~corraux_checker(input1))
    error('* corr_test2wass : incorrect input1. Please use CORR_INITIALIZE before usage.');
end
if (~corraux_checker(input2))
    error('* corr_test2wass : incorrect input2. Please use CORR_INITIALIZE before usage.');
end

%   parameters
m = input1.size(3); % number of observations from class 1
n = input2.size(3); % okay
if (input1.size(1)~=input2.size(1))
    error('* corr_test2wass : two inputs should have same dimension.');
end

%% COMPUTATION - DISTANCES
DX0 = corr_pdist(input1);
DY0 = corr_pdist(input2);
DZ0 = corraux_pdist2(input1.data, input2.data);
DXY = [DX0, DZ0;DZ0', DY0]; % concatenate for future use.

%% COMPUTATION - STATISTIC
Tmn = corraux_wassersteinD(DZ0, 2.0);

%% COMPUTATION - ITERATIONS
Tvec = zeros(1,nreps);
for i=1:nreps
    % random permutation
    id0 = randperm(m+n); % random permutation
    idx = id0(1:m);
    idy = id0((m+1):(m+n));
    
    DZ1 = DXY(idx,idy);
    Tvec(i) = corraux_wassersteinD(DZ1, 2.0);
end

%% P-VALUE
pvalue = (sum(Tvec>=Tmn)+1)/(nreps+1);
end
