function [cc, coeff] = spd_icca(input, k)

% SPD_ICCA extracts Independent Covariance Components similar to
% conventional ICA using reconstruction ICA algorithm. In return,
% transformation weights are also SPD matrices.
%
%   * USAGE
%       [pc, coeff] = spd_icca(input, k)
%   * INPUT
%       input       an object from 'spd_initialize' for (p,p,N) data.
%       k           the number of independent covariance components.
%   * OUTPUT
%       cc          (p,p,k) 3d array of independent covariance components.
%       coeff       (N,k)   matrix of coefficients.
%   * AUTHOR     Kisung You (kyou@nd.edu)
%   * HISTORY
%       0.1. [06/2019] initial implementation.

%% Preprocessing : Check Inputs
%   1. input
if (~spdaux_checker(input))
    error('* spd_iga : incorrect input. Please use SPD_INITIALIZE before usage.');
end
p = input.size(1); % dimension
N = input.size(3); % number of data points
%   2. k
k = round(k);
if ((k<=1)||(k>=p))
    error('* spd_iga : input K should be a number in (1,p).');
end

%% Compute Intrinsic Mean 
equidat = spdaux_equivariant(input); % N-by-p

%% Use 'rica'
Mdl = rica(equidat, k); 
Mdlweights = Mdl.TransformWeights; % p-by-k; transformation weights

%% wrap-up
%   1. transformation weights
cc = zeros(p,p,k);
for i=1:k
    tgt = reshape(Mdlweights(:,i),[p,p]);
    cc(:,:,i) = expm((tgt+tgt')/2);
end
%   2. coefficients
coeff = equidat*Mdlweights;


end
