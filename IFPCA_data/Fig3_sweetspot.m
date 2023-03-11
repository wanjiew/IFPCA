% This script is to show the clustering performance of IF-PCA with different 
% choices of threshold for Lung Cancer(1) data (Fig. 3 in Jin and Wang (2015)). 
%
% This script will generate a plot presenting clustering error rate vs number of 
% selected features for Lung Cancer(1) data
% 
% To run this script, it is required that the folder Data is in the path, containing 
% lungCancer.mat data set. The data set is provided at <link>


load('Data/lungCancer.mat');
Data = [lungCancer_test(1:149, 1:12533); lungCancertrain(:, 1:12533)];
%Transform the data set as p by n data set
Data = Data';
Class = [lungCancer_test(1:149, 12534); lungCancertrain(:, 12534)];
load('lungcancer_pval.mat');


n = length(Class); p = size(Data, 1);
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);

% Calculate KS statistic for Lung Cancer(1)
ind = (0:n)/n;
KS = zeros(p, 1); prob = zeros(1, p);
for j= 1:p
	pi = normcdf(Data(j,:));
	pi = sort(pi);
	KS(j) = sqrt(n)*max(max(abs(ind(1:n) - pi)), max(abs(ind(2:(n+1)) - pi)));
end
KSsort = sort(KS, 'descend');

select = 1:20:p; select = [select, 251];
errori = select*0; 
for i = 1:length(select)
	data_select = Data(KS > KSsort(select(i)), :);
	
	G = data_select'*data_select;
	% Find the leading eigenvectors
	[V, ~] = eigs(G, 1); 
	error_rep = zeros(30,1);
	for jj = 1:30
	class1a = kmeans(V, 2, 'Replicates', 30);
	% Find the error rate with the estimated labels
	t = crosstab(class1a, Class);
    error_rep(jj) = min(sum(diag(t)), n - sum(diag(t)))/n;
	end
	errori(i) = mean(error_rep);
end

figure
plot(select, errori, 'linewidth', 3)
hold on
plot([251, 251], [-0.04, max(errori)+1/n], 'r--', 'linewidth', 3);
axis([1 5600 -0.04 max(errori)+1/n ])
hold off
