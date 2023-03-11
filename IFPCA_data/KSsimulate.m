function [KSvalue] = KSsimulate(n, rep)
% This function is to simulate KS statistics with parameter $n$ 
% and the density function 
% sqrt(n) * max |Empirical CDF of x - CDF of N(mean(x), std(x))|,
% mean(x) = sum(x)/n, std(x) = sqrt(sum(x-mean(x))^2/(n-1)),
% where $x$ is n by 1 sample with i.i.d normal distribution. 
%
% Inputs: 
% n: the parameter indicating the sample size of $x$
% rep: the number of KS scores to simulate.
%
% Outputs:
% KSvalue: rep by 1 vector, containing the simulated KS scores.
%
KSvalue = zeros(rep,1); kk = (0:n)'/n;
for i = 1:rep

	x = randn(n,1); 
	z = (x - mean(x))/std(x);
	z = z/sqrt(1 - 1/n);
	pi = normcdf(z);
	pi = sort(pi);

	KSvalue(i) = max(max(abs(kk(1:n) - pi)), max(abs(kk(2:(n+1)) - pi)));
	end
KSvalue = KSvalue*sqrt(n);
end