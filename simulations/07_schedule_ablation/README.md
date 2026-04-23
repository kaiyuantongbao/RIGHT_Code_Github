# RIGHT diagnostic experiments: parallel version

Place these files under

```text
RIGHT_Code_Github/simulations/07_schedule_ablation/
```

or run them from any directory inside `RIGHT_Code_Github`. The utility file searches upward until it finds:

```text
R/solvers.R
R/data_generator.R
```

The code reuses the repository implementation of `solver_right()`, `robust_grad_MoM()`, `hard_threshold()`, data generators, and sample-wise gradient functions. It only adds diagnostic wrappers.

## Scripts

```text
31_support_localization_FD_linear_parallel.R
32_support_localization_FD_logistic_parallel.R
```

These record FD-RIGHT trajectories and produce active-set diagnostics. They output both the old cross-replication selection-frequency heatmap and additional within-trajectory summaries: support churn, moving-window union size, nuisance-window union size, and L2 error trajectory.

```text
33_gradient_complexity_pointwise_linear_parallel.R
34_gradient_complexity_pointwise_logistic_parallel.R
```

These reproduce the pointwise gradient-estimation diagnostic, but calibrate the local block count as

```text
K2 = ceil(c_K log p), c_K in {1, 2, 4, 8}.
```

This addresses the case where `K2 = ceil(log p)` is numerically too small.

```text
35_gradient_complexity_uniform_linear_parallel.R
36_gradient_complexity_uniform_logistic_parallel.R
```

These implement the fixed-iterate versus uniform-class diagnostic. For each radius, the scripts generate many candidate sparse perturbations and compute prefix maxima over `M` candidates. `M = 1` mimics fixed-iterate control; larger `M` mimics simultaneous control over a larger adaptive class. This is closer to the global-to-local complexity distinction than the pointwise radius diagnostic alone.

## Parallel execution

The code uses `parallel::mclapply()` on Unix/Linux/macOS. On Windows it falls back to sequential execution. Change `n_cores` in each script if needed.

## Output directories

Results are written under:

```text
results/07_schedule_ablation/
```

with subdirectories:

```text
E_support_localization_FD_linear/
E_support_localization_FD_logistic/
E_gradient_complexity_pointwise_linear/
E_gradient_complexity_pointwise_logistic/
E_gradient_complexity_uniform_linear/
E_gradient_complexity_uniform_logistic/
```

Each directory contains `raw/`, `summary/`, and `figures/`.

## Recommended run order

For a quick check, run:

```bash
Rscript 31_support_localization_FD_linear_parallel.R
Rscript 33_gradient_complexity_pointwise_linear_parallel.R
Rscript 35_gradient_complexity_uniform_linear_parallel.R
```

Then run the logistic counterparts after confirming the linear scripts behave as expected.

The logistic uniform diagnostic is slower because the reference gradient is computed by Monte Carlo. For a pilot run, reduce `reps`, `candidate_grid`, or `M_ref`. For final figures, increase `M_ref`.
