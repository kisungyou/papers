function output = spdaux_perturb(input, sd,mode)
% Description
%   : Given SPD matrix 'input', it adds some tangential Gaussian noise
%     of 'sd' amount/magnitude (standard deviation).
if nargin<3, mode=0; end
p = size(input,1);
if mode==0, noise = randn(p)*sd; %Gaussian noise
else, noise = rand(p)*sd; end %uniform noise

noise = (noise+noise')/2;
output = spdaux_exp(input, noise, 1.0);

end