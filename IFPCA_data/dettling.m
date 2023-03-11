%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-processing Procedure for gene expression data analysis (Dudoit et al(2000), 
%Dettling(2002)). 
%There are three steps. First, we thresholding the data by a lower bound 
%and an upper bound. Then we exclude the features with $max - min < ranb$ 
%and $max/min < ratiob$. Here $max$ and $min$ correspondingly denote the maximum 
%and minimum of this feature among the samples. At last, take base 10 logarithmic 
%transformation for the data. 
%The method is first introduced by Dudoit, and then applied by Dettling. 
%
%Function: [Data, p2] = dettling(Data, lb, ub, ranb, ratiob)
%
%Input parameters: 
%Data: p by n data matrix, where p is the number of features and n is the 
%number of samples. Each column presents the observations for a sample.  
%lb: lower bound, entries in Data smaller than this number will be 
%corrected to this bound, default to be 100
%ub: upper bound, entries in Data larger than this number will be corrected
%to this bound, default value 16000
%ranb: the bound of range. When the range of a feature is smaller than this
%bound, the feature will be excluded. The default value of it is 500
%ratiob: the bound of ratio, default value 5. When the ratio of maximum 
%value and minimum value for a feature is smaller than ratiob, the feature 
%will be excluded.
%
%Output parameters:
%p2: number of genes after cleaning
%Data: p2 by n data matrix, which is the cleaned data set. Each column of
%it presents the cleaned observations for a sample
%
%Example:
% Data = load('~/Data/SuCancer2001.txt');
% %Normalize this data with the quantiles of Leukemia Data
% Data = normquan(Data, 0.25, 0.75);
% %Clean the data with dettling's procedure and the same parameter in
% %Dettling's paper
% Data = dettling(Data);
%
%Reference:
% Dettling and Buhlmann (2002), Boosting for tumor classification with gene
%expression data
% Dudoit, Fridlyand and Speed (2000), Comparison of Discrimination Methods 
%for the Classication of Tumors Using Gene Expression Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Data, p2] = dettling(Data, lb, ub, ranb, ratiob)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%Pre-Processing%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Default values
if (nargin<2||isempty(lb))
    lb = 500;
end

if (nargin<3||isempty(ub))
    ub = 16000;
end

if (nargin<4||isempty(ranb))
    ranb = 500;
end

if (nargin<5||isempty(ratiob))
    ratiob = 5;
end

RawData = Data;
p2 = size(RawData, 1);

ind = [];
for i = 1:p2
    v = RawData(i, :);
    if sum(isnan(v)) == 0 && sum(isinf(v)) == 0
        ind = [ind, i];
    end
end
RawData = RawData(ind, :);

RawData(RawData > ub) = ub;
RawData(RawData < lb) = lb;

ind = [];
for i = 1:p2
    flag = 0;
    v = RawData(i, :);

    if max(v)/min(v) <= ratiob || range(v) <= ranb
        flag = 1;
    end
    
    if flag == 0 
        ind = [ind, i];
    end
end

Data = RawData(ind, :);
Data = log10(Data);
p2 = size(Data, 1);
 
