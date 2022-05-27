function output = corraux_perturb(input, sd, mode)

% Given a correlation matrix 'input', it adds some tangential Gaussian
% noise of 'sd' amount / magnitude.

if (nargin < 3)
    mode = 0;
end

p = size(input, 1);
if (mode==0) % case : gaussian
    noise = randn(p)*sd;
else         % case : uniform
    noise = rand(p)*sd;
end

noise  = (noise+noise')/2;
output = corraux_cov2cor(corraux_spdexp(input, noise, 1.0));


end