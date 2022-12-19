function pvalue = spd_eqtest(input1, input2, mode)

% SPD_EQTEST performs a hypothesis testing for the equality of two
% distributions given samples 'input1' and input2'. It follows the
% nonparametric method by Biswas and Ghosh (2014). 
%
%   * USAGE
%       pvalue = spd_eqtest(input1, input2)
%       pvalue = spd_eqtest(input1, input2, mode)
%   * INPUT
%       input1     an object from 'spd_initialize' for class 1
%       input2     an object from 'spd_initialize' for class 2
%       mode       type of analysis. 'extrinsic' or 'intrinsic'.
%                  Default is 'extrinsic' for faster computation.
%   * OUTPUT
%       pvalue    p-value under 'H0 : two are equally distributed'.
%
%   * AUTHOR   Kisung You (kyou@nd.edu)
%   * HISTORY
%       0.1. [03/2019] initial implementation.

%% checker and default parameters
%   checkers
if (~spdaux_checker(input1))
    error('* spd_eqtest : incorrect input1. Please use SPD_INITIALIZE before usage.');
end
if (~spdaux_checker(input2))
    error('* spd_eqtest : incorrect input2. Please use SPD_INITIALIZE before usage.');
end
if (nargin < 3)
    dmode = 'extrinsic';
else
    dmode = mode;
end
%   parameters
m = input1.size(3); % number of observations from class 1
n = input2.size(3); % okay
if (input1.size(1)~=input2.size(1))
    error('* spd_eqtest : two inputs should have same dimension.');
end
%   permutation
nreps = 1000;

%% copmute distance : DX(m,m), DY(n,n), DZ(m,n)
DX0 = spd_pdist(input1, dmode);
DY0 = spd_pdist(input2, dmode);
DZ0 = spdaux_pdist2(input1.data, input2.data, dmode);
DXY = [DX0, DZ0;DZ0', DY0]; % concatenate for future use.

Tmn = spd_eqdist_statistic(DX0,DY0,DZ0); % out test statistic

%% ASYMPTOTIC VS PERMUTATION
%% 1. PERMUTATION
Tvec = zeros(1,nreps);
for i=1:nreps
    id0 = randperm(m+n); % random permutation
    idx = id0(1:m);
    idy = id0((m+1):(m+n));
    
    DX1 = DXY(idx,idx);
    DY1 = DXY(idy,idy);
    DZ1 = DXY(idx,idy);
    Tvec(i) = spd_eqdist_statistic(DX1,DY1,DZ1);
end
pvalue = (sum(Tvec>Tmn)+1)/(nreps+1);

% 2. ASYMPTOTICS
% lbd = m/(m+n);
% S1 = spd_eqdist_computeS(DX0);
% S2 = spd_eqdist_computeS(DY0);
% sig2 = (m*S1 + n*S2)/(m+n);
% Tmnstar = ((m+n)*lbd*(1.0-lbd)/(2*sig2))*Tmn;
% pvalue = chi2cdf(Tmnstar,1,'upper');

end

%% auxiliary functions
function output = spd_eqdist_statistic(DX,DY,DXY)
m = size(DXY,1);
n = size(DXY,2);

muff = sum(sum(triu(DX,1)))/(m*(m-1)/2);
mufg = sum(DXY)/(m*n);
mugg = sum(sum(triu(DY,1)))/(n*(n-1)/2);

vec1 = [muff,mufg];
vec2 = [mufg,mugg];
output = sum((vec1-vec2).^2);
end

function output = spd_eqdist_computeS(D)
m = size(D,1);
% triplet
term1 = 0.0;
for i=1:(m-2)
    for j=(i+1):(m-1)
        for k=(j+1):m
            term1 = term1 + D(i,j)*D(i,k);
        end
    end
end
term1 = term1/((m*(m-1.0)*(m-2.0))/6.0);
% pair with sqaured
term2 = 0.0;
for i=1:(m-1)
    for j=(i+1):m
        term2 = term2 + D(i,j);
    end
end
term2  = term2/(m*(m-1)/2);
output = term1-(term2*term2);
end
