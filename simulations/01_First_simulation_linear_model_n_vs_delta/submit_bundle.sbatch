#!/bin/bash

#SBATCH --job-name=final_run_bundle
#SBATCH --output=logs/bundle_run_%A_%a.out
#SBATCH --error=logs/bundle_run_%A_%a.err

# --- [KEY] Submit 10 "bundle" jobs ---
# Total tasks = 4200. Each bundle handles ~420 tasks.
#SBATCH --array=1-10

# --- [KEY] Resource Allocation for EACH bundle job ---
# Time: (420 tasks / 20 cores) * ~30s/task = 630s = ~10.5 mins. Requesting 45 mins is safe.
#SBATCH --time=00:45:00
# Cores: Request 20 cores for internal parallelism
#SBATCH --cpus-per-task=20
# Memory: Request 4G per core
#SBATCH --mem-per-cpu=4G
# Each array task is a single process
#SBATCH --ntasks=1

#SBATCH --partition=cpu

# --- Environment and Execution ---
export WORK_DIR="/path/to/your/RIGHT_Code_Github"
cd $WORK_DIR

module load miniconda3
source activate your_r_env_name

echo "--- Starting Bundle Job ---"
echo "Job ID: $SLURM_JOB_ID, Bundle ID (Task ID): $SLURM_ARRAY_TASK_ID"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"
echo "---------------------------"

# Run the new bundle script
Rscript Rscript simulations/01_First_simulation_linear_model_n_vs_delta/02_run_hpc_bundle.R

echo "--- Bundle Job Finished ---"