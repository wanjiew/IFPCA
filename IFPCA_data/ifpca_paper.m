function [err, stats, numselect] = ifpca_paper(Data, Class, K, options) 
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

% Read the options
kmeansrep = options.kmeansrep;
Nrep = options.Nrep;
per = options.per;
pvalcut = options.pvalcut;
KSvalue = options.KS;
% End read

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize Data
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);

% Find the adjustment parameters for IF-HCT-PCA
rep = length(KSvalue);
KSvaluesort = sort(KSvalue, 'ascend');
adja = mean(KSvaluesort(1:round(per*rep))); adjb = std(KSvaluesort(1:round(per*rep)));

% Calculate KS value for each feature in the data set
KS = zeros(p, 1);
for j= 1:p
    pi = normcdf(Data(j,:)/sqrt(1 - 1/n));
    pi = sort(pi);
    kk = (0:n)'/n;
    KS(j) = max(max(abs(kk(1:n) - pi')), max(abs(kk(2:(n+1)) - pi')));
end

% Adjust the KS values with adjustment parameters from simulated KSvalue
[KSsort, ranking] = sort(KS, 'ascend');
KSadjust = (KS - mean(KSsort(1:round(per*p))))/std(KSsort(1:round(per*p)));
KSadjust = KS*adjb + adja;

%Calculate P-value with simulated KS values
[~, ind] = sort([KSvaluesort, KSadjust], 'ascend');
psort = ind((rep+1):end) - (1:p);
psort = 1 - psort/rep; 
pval = zeros(p,1); pval(ranking) = psort;

% Calculate HC functional at each data point
kk = (1:p)'/(1 + p);
HCsort = sqrt(p)*(kk - psort)./sqrt(kk);
HCsort  =   HCsort./sqrt(max(sqrt(n)*(kk - psort)./kk, 0) + 1 );
HC = zeros(p,1); HC(ranking) = HCsort;

% Decide the threshold
Ind = find(psort>pvalcut, 1, 'first');
ratio = HCsort;
ratio(1:Ind-1) = -Inf; ratio(round(p/2)+1:end)=-Inf;
L = find(ratio == max(ratio), 1, 'last');
numselect = L; 

% Record the statistics for every feature
stats.KS = KS; stats.HC = HC; stats.pval = pval; stats.ranking = ranking; 

% IF-PCA
data_select = Data(pval <= psort(L), :);
G = data_select'*data_select;
[V, ~] = eigs(G, K - 1); 
errors = zeros(Nrep, 1);
for jj=1:Nrep
    class1a = kmeans(V, K, 'replicates', kmeansrep);
    t = crosstab(class1a, Class);
	errors(jj) = errors.calculate(t);
end
err.ifpca = errors; 


errors = zeros(Nrep, 1); 
for jj = 1:Nrep
    class1a = kmeans(data_select', K, 'Replicates', kmeansrep);
    t = crosstab(class1a, Class);
	errors(jj) = errors.calculate(t);
end
err.ifkmeans = errors;

% Hiearchical clustering with complete package
d = pdist(data_select');
Z = linkage(d, 'complete');
c = cluster(Z, 'maxclust', K);
t = crosstab(c, Class);
err.ifcom = errors.calculate(t);

% IF-PCA-med

% set the adjustment parameters
adja = median(KSvaluesort(1:round(per*rep))); adjb = mad(KSvaluesort(1:round(per*rep)), 1);

% Adjust KS scores with median and MAD
KSadjust = (KS - median(KSsort(1:round(per*p))))/mad(KSsort(1:round(per*p)), 1);
KSadjust = KS*adjb + adja;

%Calculate P-value with simulated KS values
[~, ind] = sort([KSvaluesort, KSadjust], 'ascend');
psort = ind((rep+1):end) - (1:p);
psort = 1 - psort/rep; 
pval = zeros(p,1); pval(ranking) = psort;

% Calculate HC functional at each data point
kk = (1:p)'/(1 + p);
HCsort = sqrt(p)*(kk - psort)./sqrt(kk);
HCsort  =   HCsort./sqrt(max(sqrt(n)*(kk - psort)./kk, 0) + 1 );
HC = zeros(p,1); HC(ranking) = HCsort;

% Decide the threshold
Ind = find(psort>pvalcut, 1, 'first');
ratio = HCsort;
ratio(1:Ind-1) = -Inf; ratio(round(p/2)+1:end)=-Inf;
L = find(ratio == max(ratio), 1, 'last');
numselect = L; 

data_select = Data(pval <= psort(L), :);
G = data_select'*data_select;
[V, ~] = eigs(G, K - 1); 
errors = zeros(Nrep, 1);
for jj=1:Nrep
    class1a = kmeans(V, K, 'replicates', kmeansrep);
    t = crosstab(class1a, Class);
	errors(jj) = errors.calculate(t);
end
err.med = errors; 
