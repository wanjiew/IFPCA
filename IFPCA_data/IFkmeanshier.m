Nrep = 30; kmeansrep = 30;
psort = sort(stats.pval, 'ascend');
data_select = Data(stats.pval <= psort(numselect), :);
% K-means 
errors = zeros(Nrep, 1); 
for jj = 1:Nrep
    class1a = kmeans(data_select', K);
    t = crosstab(class1a, Class);
    if K == 2
        errors(jj) = min(sum(diag(t)), n - sum(diag(t)))/n;
    elseif K == 3
        errors(jj) = cluster3(t);
    elseif K == 4
        errors(jj) = cluster4(t);
    elseif K == 5
        errors(jj) = cluster5(t);
    end
end
err.ifkmeans = errors;

% Hiearchical clustering with complete package
d = pdist(data_select');
Z = linkage(d, 'complete');
c = cluster(Z, 'maxclust', K);
t = crosstab(c, Class);
if K == 2
    err.ifcom = min(sum(diag(t)), n - sum(diag(t)))/n;
elseif K == 3
    err.ifcom = cluster3(t);
elseif K == 4
    err.ifcom = cluster4(t);
elseif K == 5
    err.ifcom = cluster5(t);
end


