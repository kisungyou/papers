function [opth, cverror] = spd_smoothcv(input, x, tgth)

% SPD_SMOOTHCV performs cross-validation for optimal bandwidth selection
% across the grid of possible values.
%
%   * USAGE
%       [opth, cverror] = spd_smoothcv(input)
%       [opth, cverror] = spd_smoothcv(input, x)
%       [opth, cverror] = spd_smoothcv(input, x, tgth)
%   * INPUT
%       input       an object from 'spd_initialize' for (p,p,N) data.
%       x           1d covariate data corresponding to input.data values.
%       tgth        a vector of potential bandwidth values.
%   * OUTPUT
%       opth        an optimal bandwidth value.
%       cverror     a table of tested bandwidth values and CV errors.
%   * AUTHOR     Kisung You (kyou@nd.edu)
%   * HISTORY
%       0.1. [06/2019] initial implementation.



%% Preprocessing
if (~spdaux_checker(input))
    error('* spd_smoothcv : incorrect input. Please use SPD_INITIALIZE before usage.');
end
p = input.size(1); % dimension
N = input.size(3); % number of data points
y = spdaux_equivariant(input); % (N-by-p^2)
if (nargin < 2) || isempty(x)
    x = 1:N;
end
if (length(x)~=N)
    error('* spd_smoothcv : X should be of same length as the number of observations.');
end
if (nargin < 3)
    tgth = exp(linspace(-4,2,100));
end

%% Let's try manual checking of CV error
nh = length(tgth);
errors = zeros(1,nh);
ind = crossvalind('Kfold',x,5); % it returns a grouping variable.
for i=1:nh       
    errors(i) = spdaux_smoothcvsingle(input, x, tgth(i), ind);
end

%% finalize
idmin = find(errors==min(errors));
if (length(idmin)>1)
    idmin = idmin(1);
end
opth   = tgth(idmin); % optimal bandwidth
cverror = table(tgth, errors, 'VariableNames', {'bandwidths','errors'});
    
end
