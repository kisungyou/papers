function outcell = spdaux_indexer(label)

% corresponding to ulabel, save the index for each cluster.

ulabel = sort(unique(label));
K      = length(ulabel);

outcell = cell(K,1);
for k=1:K
    outcell{k} = find(label==ulabel(k));
end

end