%%%Normalize Data with lp quantile and up quantile. The quantiles can be
%%%changed. It is required that the data Leukemia is in current path.
%
%Inputs: 
%Data: p by n data, where p is the number of features and n is the 
%number of samples. Each column presents the observations for a sample.
%lp: The lower quantile used to normalize.
%up: The upper quantile used to normalize.
%
%Outputs: 
%Data: Normalized p by n data by lp and up quantiles.
%
%Example:
% Data = load('~/Data/SuCancer2001.txt');
% %Normalize this data with the quantiles of Leukemia Data
% Data = normquan(Data, 0.25, 0.75);
% %Clean the data with dettling's procedure and the same parameter in
% %Dettling's paper
% Data = dettling(Data);

function Data = normquan(Data, lp, up)
p = size(Data, 1); n = size(Data, 2);

dataarray = reshape(Data, n*p, 1);
Q1_data = quantile(dataarray, lp);
Q3_data = quantile(dataarray, up);
r_data = Q3_data - Q1_data;

load('Data/leukemia.mat');
p2 = size(Leukemia, 1); n2 = size(Leukemia, 2);
Leukarray = reshape(Leukemia, n2*p2, 1);
Q1_leuk = quantile(Leukarray, lp);
Q3_leuk = quantile(Leukarray, up);
r_leuk = Q3_leuk - Q1_leuk;
Data = Data/r_data*r_leuk;


