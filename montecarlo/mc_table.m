function T = mc_table(file)
%MC_TABLE  The Monte Carlo table of the post, with paired-bootstrap 95% confidence intervals.
%
%   mc_table          uses montecarlo/results/mc_estimates.mat (shipped, or rebuilt by mc_collect)
%   mc_table(file)    another file written by mc_collect
%   T = mc_table(...) also returns the numbers
%
%   For each design and object the error of the posterior median is e_r = median_r - truth in
%   replication r, and the reported ratio is
%       MSE(svk)/MSE(hom) = mean_r(e_r(svk)^2) / mean_r(e_r(hom)^2),
%   below one when the model with stochastic volatility is more precise. The table of the post shows
%   the median of the ratio over 20 objects (all columns of mc_objects.m except the two impact
%   responses to a supply shock), the two slopes beta and alpha, and the response of employment to
%   a demand shock at 4 quarters. Confidence intervals: percentile intervals from 4,000 bootstrap
%   resamples of the replications, paired across the two estimators and the 20 objects (the
%   interval of the median over objects takes the median in each resample).
mcdir = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(file), file = fullfile(mcdir, 'results', 'mc_estimates.mat'); end
L = load(file, 'est'); est = L.est;
names = est.names; z = est.truth(:)';
iuse = setdiff(1:22, [12 18]);                 % the 20 objects (impact responses to supply excluded)
icol = [1 2 16];                               % beta, alpha, employment <- demand at 4 quarters
NBOOT = 4000;
rng_state = rng;
rng(20261007, 'twister');
nd = numel(est.designs);
T = struct('designs', {est.designs}, 'labels', {est.labels}, 'names', {names}, 'R', est.R);
T.ratio = zeros(nd, 22); T.ci = zeros(nd, 22, 2); T.med20 = zeros(nd, 1); T.med20_ci = zeros(nd, 2);
for d = 1:nd
    X = est.(est.designs{d}); R = size(X.hom.med, 1);
    en = (X.svk.med - z).^2;                   % R x 22 squared errors, stochastic volatility
    ed = (X.hom.med - z).^2;                   % R x 22 squared errors, homoskedastic
    T.ratio(d,:) = mean(en, 1)./mean(ed, 1);
    bidx = randi(R, R, NBOOT);                 % the same resamples for every object and estimator
    rb = zeros(NBOOT, 22);
    for j = 1:22
        ej = en(:,j); dj = ed(:,j);
        rb(:,j) = (mean(ej(bidx), 1)./mean(dj(bidx), 1))';
    end
    T.ci(d,:,:) = reshape(quantile(rb, [0.025 0.975])', [1 22 2]);
    T.med20(d) = median(T.ratio(d, iuse));
    T.med20_ci(d,:) = quantile(median(rb(:, iuse), 2), [0.025 0.975]);
end
rng(rng_state);

fprintf('\nMSE of the posterior median, with stochastic volatility over without it (R = %d samples per row)\n', est.R);
fprintf('%-16s %-24s %-24s %-24s %-24s\n', 'Volatility', 'Median over 20 objects', 'Demand slope beta', ...
    'Supply slope alpha', 'Employment <- demand, 4q');
for d = 1:nd
    c = cell(1, 4);
    c{1} = sprintf('%.2f [%.2f, %.2f]', T.med20(d), T.med20_ci(d,:));
    for k = 1:3
        j = icol(k); c{k+1} = sprintf('%.2f [%.2f, %.2f]', T.ratio(d,j), T.ci(d,j,1), T.ci(d,j,2));
    end
    fprintf('%-16s %-24s %-24s %-24s %-24s\n', est.labels{d}, c{:});
end
fprintf('Point estimates as in the table of the post; 95%% paired-bootstrap intervals in brackets.\n');

fprintf('\nAll objects: MSE ratio [95%% interval]   (* = impact response to supply, not among the 20)\n');
fprintf('%-16s %8s', 'object', 'truth');
for d = 1:nd, fprintf('   %-22s', est.labels{d}); end
fprintf('\n');
for j = 1:22
    mk = ' '; if ~ismember(j, iuse), mk = '*'; end
    fprintf('%-15s%s %8.4f', names{j}, mk, z(j));
    for d = 1:nd, fprintf('   %.2f [%.2f, %.2f]      ', T.ratio(d,j), T.ci(d,j,1), T.ci(d,j,2)); end
    fprintf('\n');
end
for d = 1:nd
    X = est.(est.designs{d});
    fprintf('%s: chains run longer than 40,000 iterations by the chain-length rule: hom %d, svk %d of %d\n', ...
        est.labels{d}, nnz(X.hom.nsave > 30000), nnz(X.svk.nsave > 30000), est.R);
end
end
