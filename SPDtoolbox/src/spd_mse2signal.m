function mse = spd_mse2signal(input1, input2)

% SPD_MSE2SIGNAL computes mean squared error for two given signals 'input1'
% and 'input2' of same size.
%
%   * USAGE
%       mse = spd_mse2signal(input1,input2)
%   * INPUT
%       input1      an object from 'spd_initialize' for (p,p,N) data.
%       input2      an object from 'spd_initialize' for (p,p,N) data.
%   * OUTPUT
%       mse         mean squared error.
%   * AUTHOR     Kisung You (kyou@nd.edu)
%   * HISTORY
%       0.1. [07/2019] initial implementation.


%% Preprocessing : checkers
%   1. input
if ((~spdaux_checker(input1))||(~spdaux_checker(input2)))
    error('* spd_mse2signal : incorrect input. Please use SPD_INITIALIZE before usage.');
end
p = input1.size(1); % dimension
N = input1.size(3); % number of data points
if (p~=input2.size(1))
    error('* spd_mse2signal : mismatching size.');
end
if (N~=input2.size(3))
    error('* spd_mse2signal : mismatching number of copies.');
end

%% Main Computation
mse = 0; 
for i=1:N
    tgt1 = input1.data(:,:,i);
    tgt2 = input2.data(:,:,i);
    mse  = mse + spdaux_dist(tgt1,tgt2,'intrinsic')^2;
end
mse = mse/N;

end
