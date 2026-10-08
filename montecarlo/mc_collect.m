function est = mc_collect(repdir, outfile)
%MC_COLLECT  Gathers the posterior medians of the replication files into one file for mc_table.
%
%   mc_collect                     reads montecarlo/reps, writes montecarlo/results/mc_estimates.mat
%   mc_collect(repdir, outfile)    other locations
%
%   For each design (dgp1, dgp2k, dgp3) it uses replications 1..R, with R the largest number such
%   that rep_0001 ... rep_R exist for all three designs (the designs share their random numbers).
%   Saved struct est:
%     designs, labels : design names and descriptions
%     names           : names of the 22 objects (mc_objects.m)
%     truth           : 1 x 22 true values
%     R               : replications per design
%     <design>.hom.med, <design>.svk.med     : R x 22 posterior medians
%     <design>.hom.nsave, <design>.svk.nsave : R x 1 draws kept after the burn-in (chain length)
mcdir = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(repdir), repdir = fullfile(mcdir, 'reps'); end
if nargin < 2 || isempty(outfile), outfile = fullfile(mcdir, 'results', 'mc_estimates.mat'); end
addpath(fileparts(mcdir), mcdir);
L = load(fullfile(mcdir, 'calib.mat'), 'C'); C = L.C;
designs = {'dgp1', 'dgp2k', 'dgp3'};
arms = {'hom', 'svk'};
R = 0;
while all(cellfun(@(d) isfile(fullfile(repdir, d, sprintf('rep_%04d.mat', R + 1))), designs)), R = R + 1; end
assert(R > 0, 'mc_collect: no complete replications in %s', repdir);
[~, names] = mc_objects(zeros(2,0), zeros(2,17,0), zeros(2,0));
est = struct('designs', {designs}, 'labels', {{'Constant', 'Common factor', 'Shock-specific'}}, ...
    'names', {names}, 'truth', C.Ztrue, 'R', R);
for d = 1:numel(designs)
    for e = 1:2
        X = struct('med', zeros(R, 22), 'nsave', zeros(R, 1));
        for r = 1:R
            L = load(fullfile(repdir, designs{d}, sprintf('rep_%04d.mat', r)), 'S'); S = L.S;
            assert(strcmp(S.task.dgp, designs{d}) && S.task.rep == r, 'mc_collect: unexpected task in %s rep %d', designs{d}, r);
            X.med(r,:) = S.(arms{e}).Q(:,3)';
            X.nsave(r) = S.(arms{e}).nsave;
        end
        est.(designs{d}).(arms{e}) = X;
    end
end
if ~isfolder(fileparts(outfile)), mkdir(fileparts(outfile)); end
save(outfile, 'est');
fprintf('mc_collect: %d replications per design from %s saved in %s\n', R, repdir, outfile);
end
