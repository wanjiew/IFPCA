function [err, stats, numselect] = ifpca_paper_v(Data, Class, K, KSvalue, pvalcut, kmeansrep, Nrep, per) 
%Output: 
%
%s: number of features selected in IF-PCA
%err: error rate for IFPCA(err.if), kmeans(err.kmeans), 
%kmeans++(err.kmeanspp), Hiearchical clustering with 'complete' package
%(err.com)
%stats: statistics for each feature, including KS scores(stats.KS),
%P-values(stats.pval), HC functionals(stats.HC), feature ranks(stats.ranking)
%numselect: number of selected features in IF-PCA

[p, n] = size(Data);

if (nargin<4||isempty(KSvalue))
    nullsimu = true; 
else 
    nullsimu = false;
end

if (nargin<5)||isempty(pvalcut)
    pvalcut = log(p)/p;
end

if (nargin<6)||isempty(kmeansrep)
    kmeansrep = 30;
end

if (nargin<7)||isempty(Nrep)
    Nrep = 30;
end

if (nargin<8)||isempty(per)
    per = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize Data
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);

%Simulate KS values
if (nullsimu)
    rep = 100*p;
    KSvalue = zeros(rep,1); 
    for i = 1:rep
        x = randn(n,1); 
	    z = (x - mean(x))/std(x);
	    z = z/sqrt(1 - 1/n);
	    pi = normcdf(z);
	    pi = sort(pi);
        kk = (0:n)'/n;
	    KSvalue(i) = max(max(abs(kk(1:n) - pi)), max(abs(kk(2:(n+1)) - pi)));
    end
    KSvalue = KSvalue*sqrt(n);
end
KSmean = mean(KSvalue); KSstd = std(KSvalue);
if(per < 1)
    rep = length(KSvalue);
    KSvaluesort = sort(KSvalue, 'ascend');
    KSmean = mean(KSvaluesort(1:round(per*rep)));
    KSstd = std(KSvaluesort(1:round(per*rep)));
end
%Calculate KS value for each feature in the data set
KS = zeros(p, 1);
for j= 1:p
    pi = normcdf(Data(j,:)/sqrt(1 - 1/n));
    pi = sort(pi);
    kk = (0:n)'/n;
    KS(j) = max(max(abs(kk(1:n) - pi')), max(abs(kk(2:(n+1)) - pi')));
end
if(per == 1)
    KS = (KS - mean(KS))/std(KS); 
else
    KSsort = sort(KS, 'ascend');
    KS = (KS - mean(KSsort(1:round(per*p))))/std(KSsort(1:round(per*p)));
end
%Calculate P-value with simulated KS values
KSadjust = KS*KSstd + KSmean;
pval = zeros(p,1);
for i = 1:p
    pval(i) = mean(KSvalue > KSadjust(i));
end
[psort, ranking] = sort(pval, 'ascend');


% Calculate HC functional at each data point
kk = (1:p)'/(1 + p);
HCsort = sqrt(p)*(kk - psort)./sqrt(kk);
HCsort  =   HCsort./sqrt(max(sqrt(n)*(kk - psort)./kk, 0) + 1 );
HC = zeros(p,1);
HC(ranking) = HCsort;


% Classical PCA
G = Data'*Data;
errors = zeros(Nrep, 1);
for jj=1:Nrep
    [V, ~] = eigs(G, K - 1); class1a = kmeans(V, K, 'replicates', kmeansrep);
    t = crosstab(class1a, Class);
    if K == 2
        errors(jj) = min(sum(diag(t)), n - sum(diag(t)))/n;
    elseif K == 3
        errors(jj) =cluster3(t);
    elseif K == 4
        errors(jj) =cluster4(t);
    elseif K == 5
        errors(jj) =cluster5(t);
    end
end
err.PCA = errors; 

% K-means 
errors = zeros(Nrep, 1); 
for jj = 1:Nrep
    class1a = kmeans(Data', K, 'replicates', kmeansrep);
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
err.kmeans = errors;

% K-means++  
errors = zeros(Nrep, 1); 
for jj = 1:Nrep
    class1a = kmeanspp(Data, K);
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
err.kmeansplus = errors;


% Hiearchical clustering with complete package
d = pdist(Data');
Z = linkage(d, 'complete');
c = cluster(Z, 'maxclust', K);
t = crosstab(c, Class);
if K == 2
    err.com = min(sum(diag(t)), n - sum(diag(t)))/n;
elseif K == 3
    err.com = cluster3(t);
elseif K == 4
    err.com = cluster4(t);
elseif K == 5
    err.com = cluster5(t);
end


