function output = spdaux_equivariant(input)

% compute equivariant embedding; logm
p = input.size(1); % dimension
N = input.size(3); % number of data points

% stack as rows
output = zeros(N,p*p);
for i=1:N
    tgtlg = logm(input.data(:,:,i)); % equivariant embedding
    output(i,:) = reshape(tgtlg, [1, p*p]);
end
end