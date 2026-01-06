# --- 02_run_hpc_bundle.R ---
# This script runs a "bundle" of simulations in parallel on a single compute node.
# It is designed to be called by a Slurm array job.

# --- 1. Setup: Load Libraries and Functions ---
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(dplyr))

source("R/data_generator.R")
source("R/solvers.R")
source("simulations/01_First_simulation_linear_model_n_vs_delta/rate_vs_delta_linear_model.R") # Contains run_rate_vs_moment_single

# --- 2. Get Bundle Information from Slurm ---

# Get the ID of this bundle from the Slurm array task ID
bundle_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(bundle_id)) {
  bundle_id <- 1 # Default to bundle 1 for local testing
  cat("Warning: SLURM_ARRAY_TASK_ID not found. Using default bundle_id = 1.\n")
}

# Get the number of CPU cores allocated to this task by Slurm
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(num_cores)) {
  num_cores <- detectCores() - 1 # Default to all but one core for local testing
  cat(sprintf("Warning: SLURM_CPUS_PER_TASK not found. Using %d cores locally.\n", num_cores))
}

cat(sprintf("--- Bundle %d started on %d cores ---\n", bundle_id, num_cores))

# --- 3. Determine the Workload for This Bundle ---

# Load the complete parameter grid
param_grid <- readRDS("param_grid_final.rds")
total_tasks <- nrow(param_grid)

# Define how many bundles we are splitting the work into
# This MUST match the --array range in the sbatch script!
num_bundles <- 10 
bundle_size <- ceiling(total_tasks / num_bundles)

# Calculate the start and end row index for this bundle's tasks
start_index <- (bundle_id - 1) * bundle_size + 1
end_index <- min(bundle_id * bundle_size, total_tasks)

# Check if there is any work to do
if (start_index > total_tasks) {
  cat(sprintf("Bundle %d has no tasks to run. Exiting.\n", bundle_id))
  quit(save = "no")
}

# Subset the main grid to get just the tasks for this bundle
tasks_for_this_bundle <- param_grid[start_index:end_index, ]
cat(sprintf("Bundle %d will process %d tasks (from index %d to %d).\n", 
            bundle_id, nrow(tasks_for_this_bundle), start_index, end_index))

# --- 4. Execute the Bundle in Parallel ---

# Register the parallel backend
registerDoParallel(cores = num_cores)

# Define fixed parameters
P <- 600; S_STAR <- 5; THETA_STAR <- rep(0, P); THETA_STAR[1:S_STAR] <- c(5, -5, 6, -6, 7)
S <- 10; ETA <- 0.02; T_MAX <- 300; C_K_FIXED <- 1.0

# Use foreach to run the tasks in parallel
bundle_results_list <- foreach(
  i = 1:nrow(tasks_for_this_bundle),
  .combine = 'rbind',
  # --- [KEY CORRECTION] ---
  # Explicitly load all necessary packages on each parallel worker.
  .packages = c("MASS", "mvtnorm") 
) %dopar% {
  
  # It's also a good (though sometimes redundant) practice to re-source 
  # the function files inside the loop for maximum robustness.
  source("R/data_generator.R")
  source("R/solvers.R")
  source("simulations/rate_vs_delta_linear_model.R")
  
  current_params <- tasks_for_this_bundle[i, ]
  
  run_rate_vs_moment_single(
    n = current_params$n,
    p = P,
    theta_star = THETA_STAR,
    delta = current_params$delta,
    s = S,
    eta = ETA,
    T_max = T_MAX,
    c_K = C_K_FIXED,
    seed = current_params$task_id
  )
}

# Stop the parallel cluster
stopImplicitCluster()

# --- 5. Save the Aggregated Result of This Bundle ---
output_dir <- "results/raw_bundles"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
output_filename <- file.path(output_dir, sprintf("bundle_result_%d.rds", bundle_id))
saveRDS(bundle_results_list, output_filename)

cat(sprintf("--- Bundle %d finished. Aggregated result for %d tasks saved to %s ---\n",
            bundle_id, nrow(bundle_results_list), output_filename))