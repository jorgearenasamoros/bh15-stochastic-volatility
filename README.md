# Stochastic volatility in the Baumeister–Hamilton SVAR

Replication code for the post [Stochastic volatility in the Baumeister–Hamilton SVAR](https://jorgearenasamoros.github.io/blog/stochastic-volatility-baumeister-hamilton/).

The model is the bivariate labor-market SVAR of Baumeister and Hamilton (2015) with a single common stochastic-volatility factor that scales all structural shocks:

```
A y_t = B x_{t-1} + u_t,    u_t ~ N(0, e^{h_t} D),    h_t = phi h_{t-1} + sigma eta_t
```

## How to run

In MATLAB, from this folder:

```matlab
run_all
```

`run_all` runs four chains per model (with and without stochastic volatility) from dispersed starting points, prints the convergence diagnostics, computes the impulse responses and the variance and historical decompositions, saves the eight figures of the post in `figures/` and prints every number cited in the text and in the table of the slopes. It takes about six minutes. The Monte Carlo table is produced by `montecarlo/` (see below).

Figure 1, the contours of the log likelihood for the two slopes given the posterior mean of the volatility path, is drawn by `fig_loglik.m` (called from `make_figures.m`) and saved as `figures/loglik-sv.png`. Figures 2 to 8 are `volatility-factor.png`, `elasticity-posteriors.png`, `historical-decomposition.png`, `historical-decomposition-wages.png`, `unit-irfs.png`, `one-sd-employment.png` and `acf-squared-shocks.png`.

Requirements: MATLAB R2020b or later (the figures place their legends with `Layout.Tile` in a tiled layout) and the Statistics and Machine Learning Toolbox; the Monte Carlo also needs the Optimization Toolbox. Tested with R2024b.

## Files

| File | Role |
|---|---|
| `run_all.m` | Runs everything except the Monte Carlo |
| `run_chain.m` | Runs and saves one chain |
| `gibbs_sv.m` | Gibbs sampler |
| `eval_Atarget.m` | Log posterior kernel of `A`, with `B` and `D` integrated out |
| `sample_h_common.m` | Volatility path: Kim–Shephard–Chib mixture and forward-filtering backward-sampling |
| `sample_psi_sv.m` | Draws of `phi` and `sigma^2` |
| `ksc7.m` | Kim–Shephard–Chib mixture constants |
| `read_data.m`, `set_prior.m`, `setA.m`, `logP.m` | Data, prior, contemporaneous matrix, prior density of `A` |
| `mcmc_diagnostics.m` | Split R-hat and effective sample size |
| `build_draws.m` | Pools the chains, computes impulse responses and decompositions |
| `structural_ma.m`, `compute_histdecomp.m` | Structural moving average, historical decomposition |
| `make_figures.m` | Figures and numbers |
| `fig_loglik.m` | Figure 1: contours of the log likelihood with stochastic volatility |
| `myacf.m`, `mypacf.m`, `ljungbox.m` | Autocorrelations and Ljung–Box test |

## Monte Carlo

The folder `montecarlo/` reproduces the Monte Carlo table of the post: the mean squared error of the posterior median with stochastic volatility relative to that without it, from 500 simulated samples of 200 quarters for each of three volatility processes.

**The table in a second.** The posterior medians of the 3,000 estimations (500 samples, three designs, two models) are stored in `montecarlo/results/mc_estimates.mat`, so

```matlab
cd montecarlo
mc_table
```

prints the table of the post with paired-bootstrap 95% confidence intervals, followed by the ratio for each of the 22 objects.

**From scratch.**

```matlab
cd montecarlo
run_montecarlo                               % 3 designs x 500 samples, one file per sample in montecarlo/reps/
mc_collect                                   % posterior medians -> montecarlo/results/mc_estimates_rerun.mat
mc_table('results/mc_estimates_rerun.mat')
```

`run_montecarlo` distributes the samples with `parfor` when the Parallel Computing Toolbox is installed and runs them one after the other otherwise. It saves each sample as soon as it is done and skips the samples already saved, so it can be stopped and restarted; `run_montecarlo({'dgp2k'}, 1:20)` runs a subset, and `mc_collect` then gathers replications 1 to R of the designs present in `montecarlo/reps/`. `mc_collect` writes to a separate file and never overwrites the shipped `results/mc_estimates.mat` unless that path is given as its second argument; `git checkout montecarlo/results/mc_estimates.mat` restores the shipped file. On an Intel Core i7-14650HX, one sample (both models, one core) takes about 15 seconds when it runs alone and about 30 seconds when 8 samples run in parallel; estimated from these times, the 1,500 samples take about 6 hours one after the other, or roughly 1.5 to 2 hours with 8 parallel workers. Every sample runs with fixed seeds and a single computational thread (also on parallel workers); with MATLAB R2024b on that machine, re-running samples reproduces the stored estimates bit for bit. Other versions or processors can differ in the last digits, which changes individual draws of a long chain but not the results statistically.

**Design.** Each sample is drawn from

```
A0 y_t = B0 x_{t-1} + u_t,    u_it = sqrt(d0_i s_it) z_it,    z_it ~ N(0,1),
```

with eight lags and a constant, the dates 1970:Q1–2019:Q4 and, as initial conditions, the eight observed quarters before 1970:Q1. The true slopes beta = −1.03 and alpha = 0.11, the lag coefficients and constants `B0` and the common log-volatility path `hbar_t` used below come from a preliminary estimation of the stochastic-volatility model on the same data (the posterior median of the slopes and the posterior means of `B` and `h_t`). They are fixed design inputs, stored in `montecarlo/calib.mat`, and do not coincide exactly with the posterior that `run_all` produces. The true slopes lie between its posterior medians with and without stochastic volatility, (−0.99, 0.10) and (−1.12, 0.13); `B0` differs from its posterior mean of `B` with stochastic volatility by up to 0.042 in an element; and `hbar_t` is flatter than its posterior mean of `h_t` (standard deviation 0.38 against 0.59, correlation 0.99), which matters little because the amplitude of the path is set by `k` below. `d0` is the time average of the squared structural shocks `A0 y_t − B0 x_{t−1}` in the data. The three designs differ only in the variance scale `s_it`, which averages one over time, so that `d0_i` is the average variance of shock `i`:

- `dgp1`, constant: `s_it = 1`;
- `dgp2k`, common factor: `s_it = exp(k hbar_t) / mean_t exp(k hbar_t)` for both shocks, where `hbar_t` is the smoothed common log-volatility path of that preliminary estimation (mean zero); only its shape matters, because its amplitude is set by `k`;
- `dgp3`, shock-specific: `s_it = exp(k_i h_it) / mean_t exp(k_i h_it)`, where `h_it` is a smoothed stochastic-volatility estimate for each structural shock `u_it / sqrt(d0_i)` in the data, with the structural parameters fixed at the truth.

A smoothed path is flatter than the path that generated it, so taking `hbar_t` or `h_it` as the truth would make the simulated data less heteroskedastic than the actual data. The amplitudes are therefore chosen by indirect inference: a stochastic-volatility smoother with the structural parameters fixed at the truth, applied to samples simulated with amplitude `k`, gives a path whose standard deviation (median over simulated samples) equals the standard deviation of the path that the same smoother gives on the data. This yields `k = 1.61` for the common factor and `k = (2.40, 1.61)` for the demand and supply shocks. The calibration is provided as fixed input data in `montecarlo/calib.mat`; the code in this repository reads it but does not regenerate it.

The three designs use the same normal draws `z` in each replication, and each sample is estimated by both models (`hom`, homoskedastic, and `svk`, with the common volatility factor, in the code and the files) with the priors and the sampler of the main code (`mc_estimate.m`, the function form of `gibbs_sv.m`, started at the posterior mode of `A`): chains of 40,000 iterations with 10,000 of burn-in. If the effective sample size of alpha is below 800, the same chain is run to 70,000 iterations and, if it is still below 800, to the length that the observed effective sample size per draw implies for 1.3 × 800, between 100,000 and 250,000 iterations. This applies to 32 of the 3,000 chains.

The 22 objects (`mc_objects.m`) are the two slopes, the four sums of lag coefficients and the two constants of the reduced form, the unit responses of wages and employment to each shock at 0, 4 and 20 quarters, and the two variance shares at 16 quarters. The median over objects in the table uses 20 of them, because the impact responses to a supply shock are implied by those to a demand shock.

`montecarlo/calib.mat` holds the struct `C`:

| Field | Size | Content |
|---|---|---|
| `T` | 1 × 1 | sample size, 200 |
| `time` | 200 × 1 | dates of the sample, 1970.00 to 2019.75 |
| `y0` | 8 × 2 | initial conditions: observed wage and employment growth, 1968:Q1–1969:Q4 |
| `a0` | 2 × 1 | true [beta; alpha] = [−1.0256; 0.1124] |
| `A0` | 2 × 2 | true contemporaneous matrix, [−beta 1; −alpha 1] |
| `B0` | 2 × 17 | true lag coefficients and constants of the structural form |
| `Phi0` | 2 × 17 | true reduced form, `A0 \ B0` |
| `d0` | 1 × 2 | average variances of the demand and supply shocks, (0.810, 0.106) |
| `hbar` | 200 × 1 | smoothed common log-volatility path in the data, from the preliminary estimation (mean zero); shape of the `dgp2k` factor |
| `k2` | 1 × 1 | amplitude of the common factor, 1.61 |
| `s2k` | 200 × 1 | variance scale of `dgp2k`, `exp(k2 hbar) / mean(exp(k2 hbar))` |
| `h3` | 200 × 2 | smoothed log-volatility of the demand and supply shocks in the data (mean zero) |
| `k3` | 1 × 2 | amplitudes of the shock-specific paths, (2.40, 1.61) |
| `s3` | 200 × 2 | variance scales of `dgp3`, `exp(k3_i h3_i) / mean(exp(k3_i h3_i))` |
| `Ztrue` | 1 × 22 | true values of the 22 objects, `mc_objects(a0, B0, d0')` |
| `seed_z` | 1 × 1 | the normals of sample `r` use the seed `seed_z + r` |
| `seed_chain` | 1 × 1 | the chains of sample `r` use the seeds `seed_chain + 10r + 1` (no SV) and `seed_chain + 10r + 2` (SV) |

| File in `montecarlo/` | Role |
|---|---|
| `run_montecarlo.m` | Runs the Monte Carlo (resumable, `parfor` if available) |
| `mc_one_rep.m` | One sample: simulation, both estimations, chain-length rule, posterior quantiles |
| `mc_simulate.m` | Simulates a sample under each design |
| `mc_estimate.m` | Gibbs sampler, with or without stochastic volatility, as a function |
| `mc_objects.m` | The 22 objects for each posterior draw |
| `ess_geyer.m` | Effective sample size (Geyer's initial monotone sequence) |
| `post_val.m` | Log posterior of `A` with `B` and `D` integrated out, for the starting value |
| `mc_collect.m` | Gathers the posterior medians in `reps/` into `results/mc_estimates_rerun.mat` |
| `mc_table.m` | The table of the post, with bootstrap confidence intervals |
| `calib.mat` | True parameters and volatility paths |
| `results/mc_estimates.mat` | Posterior medians of the 22 objects and chain lengths, for every sample and model |

The model functions shared with the main code (`setA.m`, `eval_Atarget.m`, `logP.m`, `sample_h_common.m`, `sample_psi_sv.m`, `structural_ma.m`) are taken from the repository root.

## Data

`data/labor_data.csv`, quarterly, 1947:Q1–2026:Q2, downloaded from FRED on 3 October 2026:

- real hourly compensation, nonfarm business sector (`COMPRNFB`)
- total nonfarm payroll employment, last month of each quarter (`PAYEMS`)

The estimation sample is 1970:Q1–2019:Q4.

## Credits

`logP.m`, `setA.m`, `read_data.m`, `set_prior.m`, `montecarlo/post_val.m` and the structure of the sampler are adapted from the replication files of Baumeister, C. and J. D. Hamilton (2015), "Sign restrictions, structural vector autoregressions, and useful prior information", *Econometrica*, 83(5), 1963–1999.

## License

MIT, see `LICENSE`. The files adapted from the Baumeister and Hamilton replication code keep their attribution.
