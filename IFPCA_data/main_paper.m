% %Load the data set
Data = load('Data/brain.x.txt');
Class = load('Data/brain.y.txt');
'brain'
K = max(Class) - min(Class) + 1;
[BrainErr, ~, BrainL] = ifpca_paper(Data, Class, K);


Data = load('Data/WangBreast2005.txt');
Class = [zeros(183,1); ones(93,1)];
'breast'
K = 2;
[WangErr, ~, WangL] = ifpca_paper(Data, Class, K);

%Colon data
Data = load('Data/colon.x.txt');
Class = load('Data/colon.y.txt');
'colon'
K = 2;
[Colonerr, ~, ColonL] = ifpca_paper(Data, Class, K);

%Luekemia data set
Data = load('Data/leukemia.x.txt'); 
Class = load('Data/leukemia.y.txt'); 
'leukemia'
K = 2;
[Leukerr, ~, LeukL] = ifpca_paper(Data, Class, K);


load('Data/lungCancer.mat');
Data = [lungCancer_test(1:149, 1:12533); lungCancertrain(:, 1:12533)];
%Transform the data set as p by n data set
Data = Data';
Class = [lungCancer_test(1:149, 12534); lungCancertrain(:, 12534)];
'lung_cancer'
K = 2;
[lungerr, ~, lungL] = ifpca_paper(Data, Class, K);


Data = load('Data/BhattacherjeeLung2001.txt');
Class = [zeros(139,1); ones(64,1)];
'lung2001'
K = 2;
[Lung2001err, ~, Lung2001L] = ifpca_paper(Data, Class, K);


Data = load('Data/lymphoma.x.txt');
Class = load('Data/lymphoma.y.txt');
'lymphoma'
K = max(Class) - min(Class) + 1;
[LymphomaErr, ~, LymphomaL] = ifpca_paper(Data, Class, K);


Data = load('Data/prostate.x.txt');
Class = load('Data/prostate.y.txt');
'prostate'
K = 2;
[ProstateErr, ~, ProstateL] = ifpca_paper(Data, Class, K);


Data = load('Data/srbct.x.txt');
Class = load('Data/srbct.y.txt');
K = max(Class) - min(Class) + 1;
'srbct'
[SrbctErr, ~, SrbctL] = ifpca_paper(Data, Class, K);


Data = load('Data/SuCancer2001.txt');
Class = [zeros(83,1); ones(91,1)];
%Normalize this data with the quantiles of Leukemia Data
Data = normquan(Data, 0.25, 0.75);
%Clean the data with dettling's procedure and the same parameter in
%Dettling's paper
Data = dettling(Data, 100, 16000, 500, 5);
'sucancer'
K = max(Class) - min(Class) + 1;
[SuErr, ~, SuL] = ifpca_paper(Data, Class, K, [], [], [], [], [], [], [], 0.5);
