function [pc, coeff, var] = spd_pga(input, k)

% SPD_PGA is an implementation of Fletcher's Principal Geodesic Analysis
% algorithm, an extension of famous PCA to the manifold-valued data.
%
%   * USAGE
%       [pc, coeff, var] = spd_pga(input, k)
%   * INPUT
%       input       an object from 'spd_initialize' for (p,p,N) data.
%       k           the number of principal component matrices to be drawn.
%   * OUTPUT
%       pc          (p,p,k) 3d array of principal component matrices.
%       coeff       (N,k)   matrix of coefficients. 
%       var         (k)     a vector of projected variance.
%   * AUTHOR     Kisung You (kyou@nd.edu)
%   * HISTORY
%       0.1. [02/2019] initial implementation.

%% Preprocessing : Check Inputs
%   1. input
if (~spdaux_checker(input))
    error('* spd_pga : incorrect input. Please use SPD_INITIALIZE before usage.');
end
p = input.size(1); % dimension
N = input.size(3); % number of data points
%   2. k
k = round(k);
if ((k<=1)||(k>=p))
    error('* spd_pga : K should be a number in (1,p).');
end

%% Compute Intrinsic Mean 
[center, ~] = spd_mean(input);
stretched = zeros(N,p^2); % tangentialize
for i=1:N
    logged = real(spdaux_log(center, input.data(:,:,i)));
    logged = (logged + logged')/2;
    stretched(i,:) = logged(:); % stack as rows
end

%% extraction
S = stretched'*stretched;
[V,D] = eigs(S,k);
coeff = stretched*V;
var = diag(D);
pc = zeros(p,p,k);
for i=1:k
    y = reshape(V(:,i),[p,p]);
    z = (y+y')/2;
    pc(:,:,i) = spdaux_exp(center,z,1.0);
end
end
