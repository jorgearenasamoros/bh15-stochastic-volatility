function run_chain(SV_on, SEED, A_START, chain_name)
% RUN_CHAIN  Runs one chain of gibbs_sv.m and saves it to results/chain_<name>.mat.
% Keeps all draws of A and the SV hyperparameters, and every 10th draw of B, D and h.
ndraws = 110000;                 %#ok<NASGU> total iterations
nburn  = 10000;                  %#ok<NASGU> burn-in (proposal adapted here)
gibbs_sv;
CH.a = a_post; CH.phi = phi_post; CH.sig2 = sig2_post;
CH.B = B_post(:,:,1:10:end); CH.invd = invd_post(:,1:10:end); CH.h = h_post(:,1:10:end);
CH.acc = acceptance_ratio; CH.seed = SEED; CH.start = A_START;
outdir = fullfile(fileparts(mfilename('fullpath')),'results');
if ~isfolder(outdir), mkdir(outdir); end
save(fullfile(outdir, ['chain_' chain_name '.mat']), 'CH', '-v7.3');
fprintf('chain %s: acceptance %.3f\n', chain_name, acceptance_ratio);
end
