% This script is to generate the comparison of classical PCA clustering and IF-HCT-PCA
% for Lung Cancer(1) (Fig. 1 in Jin and Wang (2015)). 
%
% This script will generate a two panel plot. The left panel is the classical PCA
% clustering result (where the threshold t = 0.0000), and the right panel is the 
% IF-HCT-PCA clustering results (where 251 features are selected).
% 
% To run this script, it is required that the folder Data is in the path, containing 
% lungCancer.mat data set. The data set is provided at <link>

load('Data/lungCancer.mat')
Data = [lungCancer_test(1:149, 1:12533); lungCancertrain(:, 1:12533)];
Data = Data';
Class = [lungCancer_test(1:149, 12534); lungCancertrain(:, 12534)];
K = 2; [p, n] = size(Data);
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

G = Data'*Data;
[Cv, ~] = eigs(G, 1); 

data_select = Data(KS >= KSsort(251), :);
G = data_select'*data_select;
[IFv, ~] = eigs(G, 1); 

rind = randperm(n);
Class = Class(rind); Cv = Cv(rind); IFv = IFv(rind);
g1 = find(Class == 1); g2 = find(Class == 0);

figure
subplot(121)
plot(g1, Cv(g1), 'ro', g2, Cv(g2), 'g+', 'LineWidth', 2.5)
hold on
plot(1:n, 0*(1:n), 'b--', 'LineWidth', 1.5)
axis([0 n -.55 .12])
hold off
legend('ADCA', 'MPM', 'Location', 'SouthWest')
subplot(122)
plot(g1, IFv(g1), 'ro', g2, IFv(g2), 'g+', 'LineWidth', 2.5)
hold on
plot(1:n, 0*(1:n), 'b--', 'LineWidth', 1.5)
hold off
axis([0 n -.55 .12])
legend('ADCA', 'MPM', 'Location', 'SouthWest')