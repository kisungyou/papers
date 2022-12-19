function output = spd_smooth(input, h, x)

% SPD_SMOOTH implements extrinsic local regression using Gaussian kernels
% given a kernel bandwidth value.
%
%   * USAGE
%       output = spd_smooth(input)
%       output = spd_smooth(input, h)
%       output = spd_smooth(input, h, x)
%   * INPUT
%       input       an object from 'spd_initialize' for (p,p,N) data.
%       h           kernel bandwidth.
%       x           1d covariate data corresponding to input.data values.
%   * OUTPUT
%       output      smoothed object of SPD matrices from the input.
%   * AUTHOR     Kisung You (kyou@nd.edu)
%   * HISTORY
%       0.1. [06/2019] initial implementation.
%       0.2. [06/2019] modified to use function handle.



%% Preprocessing : checkers
%   1. input
if (~spdaux_checker(input))
    error('* spd_smooth : incorrect input. Please use SPD_INITIALIZE before usage.');
end
p = input.size(1); % dimension
N = input.size(3); % number of data points
%   2. h : regression smoothing parameter
if (nargin < 2)
    h = 1.0;
end

%% Preprocessing : truly, prepare
if nargin<3, x = 1:N; end
if (length(x)~=N)
    error('* spd_smooth : length of X should be same as the number of observations in ''input''.');
end

%% Now, let's use 'smootheval' function
myevalfun = spdaux_smootheval(input, h, x);

%% Let's evaluate
y3d = zeros(p,p,N);
for i=1:N
    y3d(:,:,i) = myevalfun(x(i));
end
output = spd_initialize(y3d);
    
end
