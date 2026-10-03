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

`run_all` runs four chains per model (with and without stochastic volatility) from dispersed starting points, prints the convergence diagnostics, computes the impulse responses and the variance and historical decompositions, saves the seven figures of the post in `figures/` and prints every number cited in the text and the table. It takes about six minutes.

Requirements: MATLAB R2020a or later and the Statistics and Machine Learning Toolbox. Tested with R2024b.

## Files

| File | Role |
|---|---|
| `run_all.m` | Runs everything |
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
| `myacf.m`, `mypacf.m`, `ljungbox.m` | Autocorrelations and Ljung–Box test |

## Data

`data/labor_data.csv`, quarterly, 1947:Q1–2026:Q2, downloaded from FRED on 3 October 2026:

- real hourly compensation, nonfarm business sector (`COMPRNFB`)
- total nonfarm payroll employment, last month of each quarter (`PAYEMS`)

The estimation sample is 1970:Q1–2019:Q4.

## Credits

`logP.m`, `setA.m`, `read_data.m`, `set_prior.m` and the structure of the sampler are adapted from the replication files of Baumeister, C. and J. D. Hamilton (2015), "Sign restrictions, structural vector autoregressions, and useful prior information", *Econometrica*, 83(5), 1963–1999.

## License

MIT, see `LICENSE`. The files adapted from the Baumeister and Hamilton replication code keep their attribution.
