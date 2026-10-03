function [Q,pval] = ljungbox(x, K)
% LJUNGBOX  Ljung-Box Q statistic at lag K and its p-value (chi-square_K).
% On squared residuals this is a test for ARCH / volatility clustering:
% small p-value => reject white noise => clustering present.
r = myacf(x, K);
T = numel(x);
Q = T*(T+2)*sum((r.^2)./(T-(1:K)'));
pval = 1 - chi2cdf(Q, K);
end
