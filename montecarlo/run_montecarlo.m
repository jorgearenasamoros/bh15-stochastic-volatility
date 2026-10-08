function run_montecarlo(designs, reps, use_parallel, outdir)
%RUN_MONTECARLO  Runs the Monte Carlo of the post: 500 simulated samples of 200 quarters for each
%   of three volatility designs, each sample estimated with and without stochastic volatility.
%
%   run_montecarlo                            the three designs, replications 1 to 500
%   run_montecarlo({'dgp2k'}, 1:20)           a subset
%   run_montecarlo(designs, reps, false)      sequential, even with the Parallel Computing Toolbox
%   run_montecarlo(designs, reps, [], dir)    writes the replication files to dir
%
%   Designs: 'dgp1' constant volatility, 'dgp2k' common volatility factor, 'dgp3' shock-specific
%   volatility (mc_simulate.m). Each replication (mc_one_rep.m) is saved as
%   montecarlo/reps/<design>/rep_NNNN.mat and skipped if that file already exists, so the run can be
%   interrupted and resumed. With the Parallel Computing Toolbox the replications are distributed
%   with parfor over the current parallel pool (or a pool with the default profile); otherwise they
%   run one after the other. Every replication runs with a single computational thread (set in
%   run_task, on the client or on the worker).
%
%   Run time: on an Intel Core i7-14650HX, one replication (both estimators, 40,000 iterations
%   each, one core) takes about 15 seconds when it runs alone and about 30 seconds when 8 run in
%   parallel. Estimated from these times, the 1,500 replications take about 6 hours one after the
%   other, roughly 1.5 to 2 hours with 8 parallel workers. The chain-length rule of mc_one_rep.m
%   runs 32 of the 3,000 chains longer.
%
%   Afterwards, from montecarlo/:
%     mc_collect                                   gathers the posterior medians in
%                                                  results/mc_estimates_rerun.mat
%     mc_table('results/mc_estimates_rerun.mat')   prints the table from them
if nargin < 1 || isempty(designs), designs = {'dgp1', 'dgp2k', 'dgp3'}; end
if ischar(designs), designs = {designs}; end
if nargin < 2 || isempty(reps), reps = 1:500; end
if nargin < 3 || isempty(use_parallel)
    use_parallel = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
end
mcdir = fileparts(mfilename('fullpath'));
rootdir = fileparts(mcdir);
if nargin < 4 || isempty(outdir), outdir = fullfile(mcdir, 'reps'); end
addpath(rootdir, mcdir);                       % model functions (setA, eval_Atarget, ...) live in the root
L = load(fullfile(mcdir, 'calib.mat'), 'C'); C = L.C;

% task list, interleaved by replication; finished replications are skipped
todo = struct('dgp', {}, 'rep', {});
for r = reps(:)'
    for d = 1:numel(designs)
        if ~isfile(fullfile(outdir, designs{d}, sprintf('rep_%04d.mat', r)))
            todo(end+1) = struct('dgp', designs{d}, 'rep', r); %#ok<AGROW>
        end
    end
end
for d = 1:numel(designs)
    if ~isfolder(fullfile(outdir, designs{d})), mkdir(fullfile(outdir, designs{d})); end
end
fprintf('[%s] %d replications to run (%s; replications %d to %d), output in %s\n', ...
    char(datetime('now')), numel(todo), strjoin(designs, ', '), min(reps), max(reps), outdir);
if isempty(todo), return; end

t0 = tic;
if use_parallel
    p = gcp();
    wait(parfevalOnAll(p, @addpath, 0, rootdir, mcdir));
    fprintf('parfor over %d workers\n', p.NumWorkers);
    parfor i = 1:numel(todo)
        run_task(C, todo(i), outdir);
    end
else
    for i = 1:numel(todo)
        run_task(C, todo(i), outdir);
    end
end
fprintf('[%s] done: %d replications in %.1f minutes\n', char(datetime('now')), numel(todo), toc(t0)/60);
end

function run_task(C, t, outdir)
% runs one replication with a single computational thread and saves it (through a temporary
% file, so a killed run leaves no incomplete rep_NNNN.mat)
f = fullfile(outdir, t.dgp, sprintf('rep_%04d.mat', t.rep));
if isfile(f), return; end
try
    nthr = maxNumCompThreads(1);
    restore = onCleanup(@() maxNumCompThreads(nthr));
catch                                          % not every kind of parallel worker supports it
end
t0 = tic;
S = mc_one_rep(C, t.dgp, t.rep);
tmp = fullfile(outdir, t.dgp, sprintf('tmp_%04d.mat', t.rep));
save(tmp, 'S');
movefile(tmp, f, 'f');
fprintf('[%s] %s rep %4d done in %.0f s\n', char(datetime('now')), t.dgp, t.rep, toc(t0));
end
