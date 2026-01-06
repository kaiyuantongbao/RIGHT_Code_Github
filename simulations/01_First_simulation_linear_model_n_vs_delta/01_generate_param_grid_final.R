# --- 01_generate_param_grid_final.R ---
DELTAS <- c(0.2, 0.4, 0.6, 0.8, 1.5, 5)
N_SEQUENCE <- round(300 * exp(0.2  * seq(1,25,by=4)))
N_REPLICATIONS <- 100

param_grid <- expand.grid(
  delta = DELTAS,
  n = N_SEQUENCE,
  rep_id = 1:N_REPLICATIONS
)
param_grid$task_id <- 1:nrow(param_grid)

saveRDS(param_grid, "param_grid_final.rds")
cat(sprintf("Final parameter grid created with %d total tasks.\n", nrow(param_grid)))