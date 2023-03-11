% This script is to generate the comparison of KS statistic from Lung Cancer(1) 
% data with simulated KS scores (Fig. 2 in Jin and Wang (2015)). 
%
% This script will generate a two panel plot. The left panel is the histogram of 
% KS scores from Lung Cancer(1) data, with a blue line as the density 
% function of KS scores without adjustment (theoretical null), and a red 
% line as the density function of KS scores with adjustment (empirical null).
% The right panel is the survival function from Lung Cancer(1) data (blue) and 
% empirical null (red). 
% 
% To run this script, it is required that the folder Data is in the path, containing 
% lungCancer.mat data set. The data set is provided at <link>

load('Data/lungCancer.mat')
Data = [lungCancer_test(1:149, 1:12533); lungCancertrain(:, 1:12533)];
Data = Data';
Class = [lungCancer_test(1:149, 12534); lungCancertrain(:, 12534)];
K = 2; [p, n] = size(Data);

% Normalize the data
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);

% Simulate KS statistics under theoretical null
rep = 200*12000;
KSvalue = zeros(rep, 1);
for i = 1:rep
	kk = (0:n)'/n;
	x = randn(n, 1);
	x = (x - mean(x))/std(x);
	pi = normcdf(x);
	pi = sort(pi);
	KSvalue(i) = max(max(abs(kk(1:n) - pi)), max(abs(kk(2:(n+1)) - pi)));
end
KSvalue = KSvalue;
clear kk x pi i rep;

% Calculate KS statistic for Lung Cancer(1)
ind = (0:n)/n;
KS = zeros(p, 1); prob = zeros(1, p);
for j= 1:p
	pi = normcdf(Data(j,:));
	pi = sort(pi);
	KS(j) = sqrt(n)*max(max(abs(ind(1:n) - pi)), max(abs(ind(2:(n+1)) - pi)));
end

figure
subplot(121)
KSadjust = (KS - mean(KS))/std(KS)*std(KSvalue)*1.2 + mean(KSvalue)+.02;
thrynull = (KSvalue - mean(KS))/std(KS)*std(KSvalue)*1.2 + mean(KSvalue)+.02;
empnull = KSvalue;
hist(KSadjust, 200)
[N, X] = hist(KSadjust, 200);
step = X(2) - X(1);
sep = (min(thrynull) - step/2):step:(max(thrynull)+step/2);
center = sep(1:(length(sep) - 1)) + step/2;
thry = zeros(length(sep)-1, 1);
for i = 1:length(thry)
	thry(i) = sum(thrynull >= sep(i) & thrynull < sep(i+1));
end
thry = thry/length(thrynull)*p;
plot(center, thry, 'b-', 'LineWidth', 3);
hold on

sep = (min(empnull) - step/2):step:(max(empnull)+step/2);
center = sep(1:(length(sep) - 1)) + step/2;
emp = zeros(length(sep)-1, 1);
for i = 1:length(emp)
	emp(i) = sum(empnull >= sep(i) & empnull < sep(i+1));
end
emp = emp/length(empnull)*p;
plot(center, emp, 'r-', 'LineWidth', 3);
hist(KSadjust, 200)
hold off
axis([0 0.5 0 800])
legend('Theoretical Null', 'Empirical Null', 'Location', 'NorthEast')
clear empnull thrynull;

subplot(122)
hold on;
KSadjust = (KS - mean(KS))/std(KS)*std(KSvalue) + mean(KSvalue);
t = min([KSadjust; KSvalue]):0.01:max([KSadjust; KSvalue]);
emp = 0*t; null = 0*t;
for i = 1:length(t)
	emp(i) = sum(KSadjust >= t(i))/p;
	null(i) = sum(KSvalue >= t(i))/length(KSvalue);
end
plot(t, emp, 'b-', 'LineWidth', 3);
plot(t, null, 'r-', 'LineWidth', 3);
hold off
legend('Lung Cancer', 'Simulation', 'Location', 'NorthEast')

clear;