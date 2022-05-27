function mat_pval = corr_eqtestelem(input1, input2, ntest)

% CORR_EQTESTELEM performs multiple testing for the equality of two sets of
% CORR matrices elementwise. Since exact distribution can't be obtained, we
% use permutation-based resampling test for comparing two distributions.
%
%   * USAGE
%       mat_pval = corr_eqtestelem(input1, input2)
%       mat_pval = corr_eqtestelem(input1, input2, ntest)
%       mat_pval = corr_eqtestelem(input1, input2, alpha)
%   * INPUT
%       input1     an object from 'corr_initialize' for class 1
%       input2     an object from 'corr_initialize' for class 2
%       ntest      the number of permutations
%   * OUTPUT
%       mat_pval   a matrix containing p-values for each element. Note that
%                  the recorded p-value is NOT FDR corrected.
%
%   * AUTHOR   Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [10/2021] initial implementation.



%% checker/params
if (~corraux_checker(input1))
    error('* corr_eqtestelem : incorrect input1. Please use CORR_INITIALIZE before usage.');
end
if (~corraux_checker(input2))
    error('* corr_eqtestelem : incorrect input2. Please use CORR_INITIALIZE before usage.');
end
if (nargin < 3)
    ntest = 500;
end
if (nargin < 4)
    alpha = 0.05;
end
p = input1.size(1);
m = input1.size(3); % number of observations from class 1
n = input2.size(3); % okay
if (input1.size(1)~=input2.size(1))
    error('* corr_eqtestelem : two inputs should have same dimension.');
end
dat_input1 = input1.data;
dat_input2 = input2.data;
dat_inputA = cat(3,dat_input1,dat_input2);

%% reference
mat_ref = corr_elemsingle(dat_input1, dat_input2);

%% bootstrap
%   1. aggregate bootstrap samples
mat_bot = zeros(p,p,ntest);
for i=1:ntest
    id0 = randperm(m+n);
    idx = id0(1:m);
    idy = id0((m+1):(m+n));
    mat_bot(:,:,i) = corr_elemsingle(dat_inputA(:,:,idx),dat_inputA(:,:,idy));
end
%   2. unnormalized p-value
mat_pval = zeros(p,p);
%   2-2. off-diagonal
for i=1:(p-1)
    for j=(i+1):p
        Tmn  = abs(mat_ref(i,j));
        Tvec = abs(mat_bot(i,j,:));
        Tpval = sum(Tmn <= Tvec)/ntest;
        mat_pval(i,j) = Tpval;
        mat_pval(j,i) = Tpval;
    end
end
end

%% auxiliary function
function diffmat = corr_elemsingle(dat1, dat2)
    % compute means
    [mean1,~] = corraux_mean3d(dat1);
    [mean2,~] = corraux_mean3d(dat2);
    
    % compute absolute difference
    diffmat = abs(mean1-mean2);
end
