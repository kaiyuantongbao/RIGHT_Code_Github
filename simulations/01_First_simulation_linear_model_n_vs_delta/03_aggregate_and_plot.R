# --- 03_aggregate_and_plot.R ---
# This script aggregates the final raw results downloaded from the HPC cluster,
# performs analysis, and generates publication-ready plots.

# --- 1. Setup: Load necessary libraries ---
# Ensure you have these installed on your local machine
#setwd("C:/Users/Think pad/Desktop/First_simulation_linear_model_n_vs_delta")
library(dplyr)
library(ggplot2)
library(purrr) # For map_dfr, a clean way to read and bind files

# --- 2. Define Path and Load Data ---
# IMPORTANT: This path MUST match where you downloaded the files.
# R prefers forward slashes, even on Windows.
raw_results_path <- "results/01_First_simulation_linear_model_n_vs_delta/raw_bundles"

# Get a list of all result files
result_files <- list.files(raw_results_path, pattern = "*.rds", full.names = TRUE)

# Check if files were found
if (length(result_files) == 0) {
  stop("No result files found in the specified directory. Please check the path: \n", raw_results_path)
}

cat(sprintf("Found %d result files. Aggregating now...\n", length(result_files)))

# Read all .rds files and combine them into a single data frame
final_results <- map_dfr(result_files, readRDS)

cat(sprintf("Successfully aggregated a total of %d simulation runs.\n\n", nrow(final_results)))

# --- 3. Save the Aggregated Data (Highly Recommended) ---
# This saves the combined data, so you don't have to re-run the aggregation.
# We save it one level up, in the main 'results' directory.
aggregated_output_dir <- "results/01_First_simulation_linear_model_n_vs_delta"
if (!dir.exists(aggregated_output_dir)) { dir.create(aggregated_output_dir) }
saveRDS(final_results, file.path(aggregated_output_dir, "final_aggregated_results.rds"))
cat("Aggregated results saved to 'results/final_aggregated_results.rds'\n\n")


# --- 4. Analyze and Summarize ---
# Calculate the mean L2 error for each experimental setting (n, delta)
summary_results <- final_results %>%
  group_by(n, p, delta) %>%
  summarise(
    mean_l2_error = mean(l2_error, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Create log-transformed columns for plotting
  mutate(
    log_n = log(n),
    log_error = log(mean_l2_error),
    delta_factor = as.factor(paste0("delta = ", delta))
  )

# --- 5. Plot the Final Convergence Rate Graph ---
# This plot shows the raw log-log relationship
plot_final_rate <- ggplot(summary_results, aes(x = log_n, y = log_error, color = delta_factor)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3) +
  labs(
    title = "Convergence Rate of RIGHT Estimator",
    subtitle = sprintf("p = %d, N_replications = %d", 
                       unique(summary_results$p), 
                       nrow(final_results) / nrow(summary_results)),
    x = "Log(Sample Size, n)",
    y = "Log(Mean L2 Error)",
    color = "Noise Parameter (delta)"
  ) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")

print(plot_final_rate)

# Save the plot to a high-quality file
#ggsave("final_rate_plot.png", plot_final_rate, width = 12, height = 8, dpi = 300)
#cat("Final rate plot saved to 'final_rate_plot.png'\n\n")


# --- 6. Quantitative Analysis: Estimate Slopes ---
# This part fits a linear model to the log-log data to estimate the slopes
slope_estimates <- summary_results %>%
  group_by(delta, delta_factor) %>%
  do(model = lm(log_error ~ log_n, data = .)) %>%
  summarise(
    slope = coef(model)[2],
    r_squared = summary(model)$r.squared,
    .groups = 'drop'
  )

cat("--- Estimated Slopes from lm(log_error ~ log_n) ---\n")
slope_estimates$Theoretical_rate<-c(-0.1304348, -0.2592593, -0.3548387,-0.4285714,-0.5,-0.5)
print(slope_estimates)




# --- 03_aggregate_and_plot.R (continued) ---

# --- Part 7: Plot Normalized Convergence Rates for Slope Comparison ---

# Normalize curves to start at 0 for direct slope comparison
# We group by delta_factor and subtract the first log_error value from all others.
normalized_summary_results <- summary_results %>%
  group_by(delta_factor) %>%
  mutate(
    normalized_log_error = log_error - first(log_error)
  )
normalized_summary_results$delta_minus_0.05<-as.factor(normalized_summary_results$delta-0.05)
# Create the normalized plot
plot_normalized_rate <- ggplot(normalized_summary_results, aes(x = log_n, y = normalized_log_error, color = delta_minus_0.05)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    x = "log(n)",
    y = "Normalized Log(Mean L2 Error)",
    color = "Noise Parameter (delta)"
  ) +
  theme_bw(base_size = 22) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 18)   
  )

#print(plot_normalized_rate)

# Save the normalized plot
#ggsave("final_normalized_rate_plot.png", plot_normalized_rate, width = 11, height = 4.5, dpi = 300)
#cat("Final normalized rate plot saved to 'final_normalized_rate_plot.png'\n")


#Change colors(final graphs)
# install.packages("ggplot2")
library(ggplot2)


actual_delta_levels <- c("0.15", "0.35", "0.55", "0.75", "1.45", "4.95")


user_colors <- c(
  "#4F2982", #  (Deep Indigo)
  "#3B82A1", #  (Steel Blue)
  "#00A08A", #  (Teal)
  "#6DB38B", #  (Seafoam Green)
  "#A98C63", #  (Camel)
  "#E2C16A" #  (Ginger Yellow)
)



selected_palette <- setNames(user_colors, actual_delta_levels)



plot_normalized_rate <- ggplot(
  normalized_summary_results, 
  aes(x = log_n, y = normalized_log_error, color = `delta_minus_0.05`) # 建议用反引号包裹非标准变量名
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  

  labs(
    x = "Log(Sample Size, n)",
    y = "Normalized Log(Mean L2 Error)",
    color = "Noise Parameter" 
  ) +
  
 
  scale_color_manual(
    name = "Noise Parameter (Delta)", 
    values = selected_palette
  ) +
  theme_bw(base_size = 22) + 
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    
    
    legend.position = c(0.003, 0.006),
    legend.justification = c(0, 0),
    
    
    legend.background = element_rect(fill = "white", color="black",linewidth = 0.5),
    legend.key = element_rect(fill = "transparent"),
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 16),
    
   
    panel.grid.minor = element_blank()
  ) +
  
 
  guides(color = guide_legend(ncol = 2))


print(plot_normalized_rate)

ggsave(file.path(aggregated_output_dir,"final_normalized_rate_plot.pdf"), plot_normalized_rate,width = 8.2, height = 4.5, dpi = 300) 
ggsave(file.path(aggregated_output_dir,"final_normalized_rate_plot.pdf"), plot_normalized_rate, width = 8.2, height = 4.5, dpi = 300)
