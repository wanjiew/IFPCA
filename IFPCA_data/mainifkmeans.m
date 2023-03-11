% Brain 
Data = load('Data/brain.x.txt');
Class = load('Data/brain.y.txt');
n = size(Data,2);
K = max(Class)-min(Class)+1;
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);
load('brainresults.mat')
IFkmeanshier;
save('brainresults.mat', 'err', 'stats', 'numselect');
[err.ifcom, mean(err.ifkmeans), std(err.ifkmeans)]


% Breast
Data = load('Data/WangBreast2005.txt');
Class = [zeros(183,1); ones(93,1)];
n = size(Data,2);
K = max(Class)-min(Class)+1;
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);
load('breastresults.mat')
IFkmeanshier;
save('breastresults.mat', 'err', 'stats', 'numselect');
[err.ifcom, mean(err.ifkmeans), std(err.ifkmeans)]

%Colon 
Data = load('Data/colon.x.txt');
Class = load('Data/colon.y.txt');
n = size(Data,2);
K = max(Class)-min(Class)+1;
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);
load('colonresults.mat')
IFkmeanshier;
save('colonresults.mat', 'err', 'stats', 'numselect');
[err.ifcom, mean(err.ifkmeans), std(err.ifkmeans)]


% Leukemia 
Data = load('Data/leukemia.x.txt'); 
Class = load('Data/leukemia.y.txt'); 
n = size(Data,2);
K = max(Class)-min(Class)+1;
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);
load('leukemiaresults.mat')
IFkmeanshier;
save('leukemiaresults.mat', 'err', 'stats', 'numselect');
[err.ifcom, mean(err.ifkmeans), std(err.ifkmeans)]


% Lung1
load('Data/lungCancer.mat');
Data = [lungCancer_test(1:149, 1:12533); lungCancertrain(:, 1:12533)];
Data = Data';
Class = [lungCancer_test(1:149, 12534); lungCancertrain(:, 12534)];
n = size(Data,2);
K = max(Class)-min(Class)+1;
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);
load('lung1results.mat')
IFkmeanshier;
save('lung1results.mat', 'err', 'stats', 'numselect');
[err.ifcom, mean(err.ifkmeans), std(err.ifkmeans)]

% Lung2
Data = load('Data/BhattacherjeeLung2001.txt'); n = size(Data,2);
Class = [zeros(139,1); ones(64,1)];
K = max(Class)-min(Class)+1;
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);
load('lung2results.mat')
IFkmeanshier;
save('lung2results.mat', 'err', 'stats', 'numselect');
[err.ifcom, mean(err.ifkmeans), std(err.ifkmeans)]


%Lymphoma
Data = load('Data/lymphoma.x.txt');
Class = load('Data/lymphoma.y.txt');
n = size(Data, 2);
K = max(Class)-min(Class)+1;
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);
load('lymphomaresults.mat')
IFkmeanshier;
save('lymphomaresults.mat', 'err', 'stats', 'numselect');
[err.ifcom, mean(err.ifkmeans), std(err.ifkmeans)]


%Prostate
Data = load('Data/prostate.x.txt');
Class = load('Data/prostate.y.txt');
n = size(Data,2);
K = max(Class)-min(Class)+1;
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);
load('prostateresults.mat')
IFkmeanshier;
save('prostateresults.mat', 'err', 'stats', 'numselect');
[err.ifcom, mean(err.ifkmeans), std(err.ifkmeans)]


%SRBCT
Data = load('Data/srbct.x.txt');
Class = load('Data/srbct.y.txt');
n = size(Data,2);
K = max(Class)-min(Class)+1;
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);
load('srbctresults.mat')
IFkmeanshier;
save('srbctresults.mat', 'err', 'stats', 'numselect');
[err.ifcom, mean(err.ifkmeans), std(err.ifkmeans)]


% Su-Cancer
Data = load('Data/SuCancer.txt');
Class = load('Data/SuCancer.y.txt');n = size(Data,2);
K = max(Class)-min(Class)+1;
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);
load('suresults.mat')
IFkmeanshier;
save('suresults.mat', 'err', 'stats', 'numselect');
[err.ifcom, mean(err.ifkmeans), std(err.ifkmeans)]
