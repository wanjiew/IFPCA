%Load data 
load('lungCancer.mat');
Data = [lungCancer_test(1:149, 1:12533); lungCancertrain(:, 1:12533)];
Data = Data';
Class = [lungCancer_test(1:149, 12534); lungCancertrain(:, 12534)];
[p, n] = size(Data);

% M-fold CV

%%%%% Clip Threshold %%%%%%%%%%%%%
rng(5) %Set the seed for the permutation
rep = 20; err = zeros(rep, 1); 
for kk = 1:rep
	idx = randperm(n);
    err(kk) = 0;
    for ii = 1:M
		testidx = idx > round(n/M)*(ii-1) & idx <= min(round(ii*n/M), n);
		trainidx = ~testidx;
		test = Data(:, testidx); ytest = Class(testidx);
		train = Data(:, trainidx); ytrain = Class(trainidx);
		[wt, stats] = HCclassification(train, ytrain, 'clip');
		[label, score] = HCclassification_fit(wt, stats.xbar, stats.s, test);
		err(kk) = err(kk) + sum(label ~= ytest);		
    end
end
result = ['Result for ', num2str(M), '-fold CV with clip threshold is ', num2str(mean(err))];
disp(result)

%%%%% Soft Threshold %%%%%%%%%%%%%
rng(5) %Set the seed for the permutation
rep = 20; err = zeros(rep, 1); 
for kk = 1:rep
	idx = randperm(n);
    err(kk) = 0;
    for ii = 1:M
		testidx = idx > round(n/M)*(ii-1) & idx <= min(round(ii*n/M), n);
		trainidx = ~testidx;
		test = Data(:, testidx); ytest = Class(testidx);
		train = Data(:, trainidx); ytrain = Class(trainidx);
		[wt, stats] = HCclassification(train, ytrain, 'soft');
		[label, score] = HCclassification_fit(wt, stats.xbar, stats.s, test);
		err(kk) = err(kk) + sum(label ~= ytest);		
    end
end
result = ['Result for ', num2str(M), '-fold CV with soft threshold is ', num2str(mean(err))];
disp(result)


%%%%% Hard Threshold %%%%%%%%%%%%%
rng(5) %Set the seed for the permutation
rep = 20; err = zeros(rep, 1); 
for kk = 1:rep
	idx = randperm(n);
    err(kk) = 0;
    for ii = 1:M
		testidx = idx > round(n/M)*(ii-1) & idx <= min(round(ii*n/M), n);
		trainidx = ~testidx;
		test = Data(:, testidx); ytest = Class(testidx);
		train = Data(:, trainidx); ytrain = Class(trainidx);
		[wt, stats] = HCclassification(train, ytrain, 'hard');
		[label, score] = HCclassification_fit(wt, stats.xbar, stats.s, test);
		err(kk) = err(kk) + sum(label ~= ytest);		
    end
end
result = ['Result for ', num2str(M), '-fold CV with hard threshold is ', num2str(mean(err))];
disp(result)
