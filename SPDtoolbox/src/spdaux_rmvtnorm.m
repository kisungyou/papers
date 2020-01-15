function samples = spdaux_rmvtnorm(n, mu, sigma)

p = length(mu);
samples = zeros(n,p);
for i=1:n
    samples(i,:) = mvnrnd(mu, sigma);
end
end