function S = mc_one_rep(C, design, rep)
%MC_ONE_REP  One replication of the Monte Carlo: simulates a sample and estimates it twice.
%
%   S = mc_one_rep(C, design, rep)
%
%   C      : calibration (montecarlo/calib.mat)
%   design : 'dgp1' (constant volatility), 'dgp2k' (common factor), 'dgp3' (shock-specific)
%   rep    : replication number (sets the seeds)
%
%   Estimators, with the same priors and sampler (mc_estimate.m):
%     'hom' : homoskedastic BH15 model,              seed C.seed_chain + 10*rep + 1;
%     'svk' : BH15 with a common volatility factor,  seed C.seed_chain + 10*rep + 2.
%   Chains: 40,000 iterations, 10,000 of burn-in, B and D kept every 5th draw. Fixed rule for the
%   chain length, applied to each estimator and based only on the chain itself: if the effective
%   sample size of alpha is below 800, the same chain (same seed and burn-in, so its first
%   iterations are identical) is run to 70,000 iterations; if it is still below 800, it is run to
%       NB + 1.3*800/(ESS per draw), rounded up to the next 10,000 and kept within [100,000, 250,000].
%
%   Output S: task (dgp, rep), cfg, s (variance scales of the sample) and, for each estimator,
%     Q    : 22 x 5 posterior quantiles (2.5, 16, 50, 84, 97.5%) of the objects of mc_objects.m
%            (beta and alpha from all draws, the others from the thinned draws)
%     mean : 22 x 1 posterior means
%     ess_ab, ess_alpha_first, nsave, NIT, retried, extended, acc_post, a_start, seed, secs
%     and, for 'svk', hq (T x 5 quantiles of h_t), phi and sig2 (mean, 2.5, 50, 97.5%).
cfg = struct('NIT', 40000, 'NB', 10000, 'thin', 5, 'NIT_retry', 70000, 'ess_target', 800, ...
    'NIT_ext', [100000 250000], 'qv', [0.025 0.16 0.5 0.84 0.975], 'tstart', 9, 'tend', 208);
[yall, s] = mc_simulate(C, design, rep);
S.task = struct('dgp', design, 'rep', rep);
S.cfg = cfg;
S.s = single(s);
arms = {'hom', 'svk'};
for e = 1:2
    SV_on = e == 2;
    opt = struct('NIT', cfg.NIT, 'NB', cfg.NB, 'thin', cfg.thin, 'seed', C.seed_chain + 10*rep + 1 + SV_on);
    t0 = tic;
    out = mc_estimate(yall, cfg.tstart, cfg.tend, SV_on, opt);
    ess = ess_geyer(out.a(2,:)');
    E = struct('retried', false, 'extended', false, 'ess_alpha_first', ess);
    if ess < cfg.ess_target                       % first extension: 70,000 iterations
        opt.NIT = cfg.NIT_retry;
        out = mc_estimate(yall, cfg.tstart, cfg.tend, SV_on, opt);
        ess = ess_geyer(out.a(2,:)');
        E.retried = true;
        if ess < cfg.ess_target                   % second extension, from the observed ESS per draw
            rate = ess/size(out.a,2);
            opt.NIT = cfg.NB + 1e4*ceil(1.3*cfg.ess_target/rate/1e4);
            opt.NIT = min(max(opt.NIT, cfg.NIT_ext(1)), cfg.NIT_ext(2));
            out = mc_estimate(yall, cfg.tstart, cfg.tend, SV_on, opt);
            E.extended = true;
        end
    end
    % ---- objects on the thinned draws (A paired with B and D) ----
    ith = out.thin:out.thin:size(out.a,2);
    Z = mc_objects(out.a(:,ith), out.B, 1./out.invd);
    Q = zeros(22, numel(cfg.qv));
    Q(1:2,:)  = quantile(out.a', cfg.qv)';        % beta, alpha: all draws
    Q(3:22,:) = quantile(Z(:,3:22), cfg.qv)';
    E.Q = Q;
    E.mean = [mean(out.a,2); mean(Z(:,3:22),1)'];
    E.ess_ab = ess_geyer(out.a');
    E.nsave = size(out.a,2);  E.NIT = opt.NIT;
    E.acc_post = out.acc_post;  E.a_start = out.a_start;  E.seed = opt.seed;
    if SV_on
        E.hq = single(quantile(out.h', cfg.qv)');
        E.phi = [mean(out.phi) quantile(out.phi, [0.025 0.5 0.975])];
        E.sig2 = [mean(out.sig2) quantile(out.sig2, [0.025 0.5 0.975])];
    end
    E.secs = toc(t0);
    S.(arms{e}) = E;
end
end
