% RUN_ALL  Replicates the figures, the table and the numbers of the post
% "Stochastic volatility in the Baumeister-Hamilton SVAR"
% (https://jorgearenasamoros.github.io/blog/stochastic-volatility-baumeister-hamilton/).
%
% 1. four chains per model (SV and homoskedastic) from dispersed starting points
% 2. convergence diagnostics (split R-hat, effective sample size)
% 3. impulse responses, variance and historical decompositions on the pooled draws
% 4. figures (figures/) and numbers (printed)

clear; clc;
starts = {[-0.6;0.6], [-0.2;0.05], [-3;0.03], [-1;1.5]};     % [beta; alpha]
for SV_on = [1 0]
    for c = 1:4
        if SV_on, name = sprintf('sv_%d',c); else, name = sprintf('nosv_%d',c); end
        run_chain(SV_on, c, starts{c}, name);
    end
end
mcmc_diagnostics;
build_draws;
make_figures;
