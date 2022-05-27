function mse = corr_mse2signal(input1, input2)

% CORR_MSE2SIGNAL computes mean squared error for two given signals 'input1'
% and 'input2' of same size.
%
%   * USAGE
%       mse = corr_mse2signal(input1,input2)
%
%   * INPUT
%       input1      an object from 'corr_initialize' for (p,p,N) data.
%       input2      an object from 'corr_initialize' for (p,p,N) data.
%
%   * OUTPUT
%       mse         mean squared error.
%
%   * AUTHOR     Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.


%% Preprocessing : checkers
%   1. input
if ((~corraux_checker(input1))||(~corraux_checker(input2)))
    error('* corr_mse2signal : incorrect input. Please use CORR_INITIALIZE before usage.');
end
p = input1.size(1); % dimension
N = input1.size(3); % number of data points
if (p~=input2.size(1))
    error('* corr_mse2signal : mismatching size.');
end
if (N~=input2.size(3))
    error('* corr_mse2signal : mismatching number of copies.');
end

%% Main Computation
mse = 0; 
for i=1:N
    tgt1 = input1.data(:,:,i);
    tgt2 = input2.data(:,:,i);
    mse  = mse + corraux_dist(tgt1,tgt2)^2;
end
mse = mse/N;

end