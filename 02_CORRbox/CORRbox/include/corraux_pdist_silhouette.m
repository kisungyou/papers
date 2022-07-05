function score = corraux_pdist_silhouette(distmat, label)

label  = round(label);
ulabel = sort(unique(label));
K      = length(ulabel);

n       = size(distmat, 1);
if (n~=K)
    error("not matching.");
end
vec_a   = zeros(n,1); % a(i)
vec_b   = zeros(n,1); % b(i)
indexer = corraux_indexer(label); % index per cluster

for i=1:n
    % a(i) : average within cluster
    %        singleton index should be taken care.
    idx1     = indexer{(ulabel==label(i))};
    if (length(idx1)==1)
        vec_a(i) = 0;
    else
        vec_a(i) = mean(distmat(i,setdiff(idx1,i)));
    end
    % b(i) : weird.. but gotta do.
    otherlabels = setdiff(ulabel, label(i)); % non-overlapping classes
    mindistance = zeros(K-1,1);
    for k=1:(K-1)
        idx2 = indexer{(ulabel==otherlabels(k))};
        mindistance(k) = mean(distmat(i,idx2));
    end
    vec_b(i) = min(mindistance);
end

% final computation
vec_s = zeros(n,1);
for i=1:n
    term_top = vec_b(i)-vec_a(i);
    term_bot = max(vec_a(i),vec_b(i));
    vec_s(i) = term_top/term_bot;
end
score = mean(vec_s);

end
