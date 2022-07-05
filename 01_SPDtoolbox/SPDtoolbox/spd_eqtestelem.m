function mat_pval = spd_eqtestelem(input1, input2, ntest, alpha)

% SPD_EQTESTELEM performs multiple testing for the equality of two sets of
% SPD matrices elementwise. Since exact distribution can't be obtained, we
% use permutation-based resampling test for comparing two distributions.
%
%   * USAGE
%       mat_pval = spd_eqtestelem(input1, input2)
%       mat_pval = spd_eqtestelem(input1, input2, ntest)
%       mat_pval = spd_eqtestelem(input1, input2, alpha)
%   * INPUT
%       input1     an object from 'spd_initialize' for class 1
%       input2     an object from 'spd_initialize' for class 2
%       ntest      the number of permutations
%       alpha      desired level of significance for FDR.
%   * OUTPUT
%       mat_pval   a matrix containing p-values for each element. Note that
%                  the recorded p-value is NOT FDR corrected.
%
%   * AUTHOR   Kisung You (kyou@nd.edu)
%   * HISTORY
%       0.1. [06/2019] initial implementation.


%% checker and default parameters
%   checkers
if (~spdaux_checker(input1))
    error('* spd_eqtestelem : incorrect input1. Please use SPD_INITIALIZE before usage.');
end
if (~spdaux_checker(input2))
    error('* spd_eqtestelem : incorrect input2. Please use SPD_INITIALIZE before usage.');
end
if (nargin < 3)
    ntest = 2000;
end
if (nargin < 4)
    alpha = 0.05;
end

%   parameters
p = input1.size(1);
m = input1.size(3); % number of observations from class 1
n = input2.size(3); % okay
if (input1.size(1)~=input2.size(1))
    error('* spd_eqtestelem : two inputs should have same dimension.');
end

% separate data
dat_input1 = input1.data;
dat_input2 = input2.data;
dat_inputA = cat(3,dat_input1,dat_input2);


%% reference
mat_ref = spd_elemsingle(dat_input1, dat_input2);

%% bootstrap
%   1. aggregate bootstrap samples
mat_bot = zeros(p,p,ntest);
for i=1:ntest
    id0 = randperm(m+n);
    idx = id0(1:m);
    idy = id0((m+1):(m+n));
    mat_bot(:,:,i) = spd_elemsingle(dat_inputA(:,:,idx),dat_inputA(:,:,idy));
end
%   2. unnormalized p-value
mat_pval = zeros(p,p);
%   2-1. diagonal
for j=1:p
    Tmn  = abs(mat_ref(j,j));
    Tvec = abs(mat_bot(j,j,:));
    mat_pval(j,j) = sum(Tmn <= Tvec)/ntest;
end
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

%% auxiliary functions
function diffmat = spd_elemsingle(dat1, dat2)
    p = size(dat1,1);
    m = size(dat1,3);
    n = size(dat2,3);
    
    % direct definition
    stdat1 = struct(); stdat1.name = 'spd'; stdat1.size = [p,p,m]; 
    stdat2 = struct(); stdat2.name = 'spd'; stdat2.size = [p,p,n];
    stdat1.data = dat1;
    stdat2.data = dat2;
    
    % compute means
    [mean1,~] = spd_mean(stdat1);
    [mean2,~] = spd_mean(stdat2);
    
    % compute absolute difference
    diffmat = abs(mean1-mean2);
end




