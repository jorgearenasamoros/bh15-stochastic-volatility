function mcmc_diagnostics()
% MCMC_DIAGNOSTICS  Split R-hat and effective sample size of beta and alpha,
% pooling the four chains of each model.
rdir = fullfile(fileparts(mfilename('fullpath')),'results');
nm = {'beta','alpha'};
for mdl = {'sv','nosv'}
    X = cell(1,4);
    for c = 1:4
        L = load(fullfile(rdir, sprintf('chain_%s_%d.mat', mdl{1}, c)), 'CH'); X{c} = L.CH.a;
    end
    for j = 1:2
        ch = cellfun(@(x) x(j,:)', X, 'UniformOutput', false);
        ess = sum(cellfun(@ess_one, ch));
        fprintf('%-5s %-6s split R-hat %.3f   effective sample size %6.0f\n', mdl{1}, nm{j}, split_rhat(ch), ess);
    end
end
end

function r = split_rhat(ch)
halves = {};
for c = 1:numel(ch)
    x = ch{c}; h = floor(numel(x)/2);
    halves{end+1} = x(1:h); halves{end+1} = x(h+1:2*h); %#ok<AGROW>
end
n = min(cellfun(@numel, halves)); Y = cell2mat(cellfun(@(x) x(1:n), halves, 'UniformOutput', false));
W = mean(var(Y,0,1)); B = n*var(mean(Y,1),0);
r = sqrt(((n-1)/n*W + B/n)/W);
end

function e = ess_one(x)
% autocorrelations by FFT, summed up to the first lag below 0.05
x = x - mean(x); n = numel(x);
f = fft(x, 2*n); ac = real(ifft(f.*conj(f))); ac = ac(1:n)/ac(1);
k = find(ac(2:end) < 0.05, 1);
if isempty(k), k = n; end
e = n/(1 + 2*sum(ac(2:k)));
end
