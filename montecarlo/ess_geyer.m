function e = ess_geyer(x)
%ESS_GEYER  Effective sample size of a chain (column vector) or of each column of a matrix.
%   Autocorrelations by FFT and Geyer's (1992) initial monotone sequence estimator:
%   Gamma_m = rho_{2m} + rho_{2m+1} are summed while positive, with monotonicity enforced;
%   ESS = N / tau, tau = -1 + 2*sum_m Gamma_m.
if isrow(x), x = x(:); end
[N, K] = size(x);
e = nan(1,K);
nf = 2^nextpow2(2*N);
for j = 1:K
    y = x(:,j) - mean(x(:,j));
    if all(y == 0) || any(~isfinite(y)), continue; end
    f  = fft(y, nf);
    ac = real(ifft(f.*conj(f)));
    ac = ac(1:N)/ac(1);
    m  = floor(N/2);
    G  = ac(1:2:2*m-1) + ac(2:2:2*m);
    ix = find(G <= 0, 1);
    if ~isempty(ix), G = G(1:ix-1); end
    G  = cummin(G);
    tau = max(-1 + 2*sum(G), 1/N);
    e(j) = N/tau;
end
end
