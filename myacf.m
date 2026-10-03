function r = myacf(x, K)
% MYACF  Sample autocorrelation function r_1..r_K (no toolbox needed).
x = x(:) - mean(x);
T = numel(x);
c0 = sum(x.^2)/T;
r = zeros(K,1);
for k = 1:K
    r(k) = (sum(x(k+1:T).*x(1:T-k))/T) / c0;
end
end
