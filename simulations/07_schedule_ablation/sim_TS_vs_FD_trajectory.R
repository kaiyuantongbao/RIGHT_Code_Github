rm(list = ls())

library(here)
library(dplyr)
library(ggplot2)
library(doParallel)
library(foreach)
library(grid)

# Use the seed-enabled two-stage schedule wrappers.
# This file assumes utils_schedule_right.R is the updated version with q-split and seed control.
source(here("simulations", "07_schedule_ablation", "utils_schedule_right.R"))

out_dir <- "results/07_schedule_ablation/TS_vs_FD_trajectory"
dir.create(file.path(out_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Experiment configuration
# ==============================================================================

cfg <- list(
  n = 800,
  p = 800,
  s_star = 5,
  s = 10,
  eta = 0.02,
  
  # Fixed-fraction TS split: n2 = floor(q * n), n1 = n - n2.
  # For this trajectory figure we vary T1 and keep q,T2 fixed.
  # Later tuning scripts will use q_grid = c(0.1, 0.2, 0.25, 0.5).
  q_grid = c(0.5),
  T1_grid = c( 150, 175,200),
  T2_grid = c(8),
  m_min = 10,
  reps = 100,
  seed_base = 20260417L,
  
  df_X = 2.5,
  scale_X = 1,
  df_eps = 1.5,
  theta_magnitude = 5
)

theta_star <- make_theta_star(cfg$p, cfg$s_star, magnitude = cfg$theta_magnitude)

schedule_grid <- expand.grid(
  q = cfg$q_grid,
  T1 = cfg$T1_grid,
  T2 = cfg$T2_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  arrange(q, T1, T2) %>%
  mutate(
    config_id = row_number(),
    Config = sprintf("TS (q=%.2f, T1=%d, T2=%d)", q, T1, T2)
  )

write.csv(
  schedule_grid,
  file.path(out_dir, "summary", "trajectory_schedule_grid.csv"),
  row.names = FALSE
)

num_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# ==============================================================================
# 2. Parallel Monte Carlo loop
# ==============================================================================

raw_results <- foreach(
  rep_id = 1:cfg$reps,
  .combine = dplyr::bind_rows,
  .export = c("cfg", "theta_star", "schedule_grid"),
  .packages = c("dplyr", "MASS", "mvtnorm", "here")
) %dopar% {
  source(here("simulations", "07_schedule_ablation", "utils_schedule_right.R"))
  
  seed_data <- cfg$seed_base + rep_id
  seed_algo_ts <- cfg$seed_base + 100000L + rep_id
  
  # Same generated data are used by all methods in this replication.
  dat <- gen_one_linear_dataset(
    seed = seed_data,
    n = cfg$n,
    p = cfg$p,
    theta_star = theta_star,
    df_X = cfg$df_X,
    Sigma_X = diag(cfg$p),
    scale_X = cfg$scale_X,
    df_eps = cfg$df_eps
  )
  
  X <- dat$X
  y <- dat$y
  
  # TS configurations.
  # The same seed_algo_ts is intentionally used across the TS configs in the
  # same replication. When q is the same, this gives the same data split and
  # isolates the effect of T1. If q differs, it still uses the same base
  # permutation before taking different split sizes.
  ans_ts <- lapply(seq_len(nrow(schedule_grid)), function(ii) {
    sc <- schedule_grid[ii, ]
    
    out <- run_ts_right_linear(
      X = X,
      y = y,
      theta_star = theta_star,
      s = cfg$s,
      eta = cfg$eta,
      q = sc$q,
      T1 = sc$T1,
      T2 = sc$T2,
      m_min = cfg$m_min,
      record_trace = TRUE,
      record_l2 = TRUE,
      record_initial = TRUE,
      seed = seed_algo_ts
    )
    
    out$Config <- sc$Config
    out$arm_type <- "TS"
    out$config_id <- sc$config_id
    out$seed_data <- seed_data
    out$seed_algo <- seed_algo_ts
    out
  })
  
  df_ts <- bind_rows(ans_ts)
  
  # FD baseline. It is run long enough to cover both the raw TS iteration axis
  # and the budget-matched FD-equivalent axis.
  max_budget <- max(df_ts$T_fd_budget, na.rm = TRUE)
  max_ts_iter <- max(schedule_grid$T1 + schedule_grid$T2)
  T_fd_run <- ceiling(max(max_ts_iter, max_budget) * 1.1)
  
  out_fd <- run_fd_right_linear(
    X = X,
    y = y,
    theta_star = theta_star,
    s = cfg$s,
    eta = cfg$eta,
    T_fd = T_fd_run,
    m_min = cfg$m_min,
    record_trace = TRUE,
    record_l2 = TRUE,
    record_initial = TRUE
  )
  
  out_fd$Config <- "FD (Baseline)"
  out_fd$arm_type <- "FD"
  out_fd$config_id <- 0L
  out_fd$seed_data <- seed_data
  out_fd$seed_algo <- NA_integer_
  
  bind_rows(df_ts, out_fd) %>%
    mutate(rep_id = rep_id)
}

parallel::stopCluster(cl)

saveRDS(
  raw_results,
  file.path(out_dir, "raw", "raw_TS_vs_FD_trajectory.rds")
)

# ==============================================================================
# 3. Summary and visualization
# ==============================================================================

trace_summary <- raw_results %>%
  filter((arm_type == "FD") | (arm_type == "TS" & eligible == TRUE)) %>%
  group_by(Config, arm_type, t) %>%
  summarise(
    median_l2 = median(l2_error, na.rm = TRUE),
    q25_l2 = quantile(l2_error, 0.25, na.rm = TRUE),
    q75_l2 = quantile(l2_error, 0.75, na.rm = TRUE),
    n_rep = dplyr::n_distinct(rep_id),
    .groups = "drop"
  )

budget_lines <- raw_results %>%
  filter(arm_type == "TS", eligible == TRUE) %>%
  group_by(Config) %>%
  summarise(
    T_budget = first(T_fd_budget),
    fd_equiv_final = max(fd_equiv_iter, na.rm = TRUE),
    q = first(q),
    T1 = first(T1),
    T2 = first(T2),
    n1 = first(n1),
    n2 = first(n2),
    b2 = first(b2),
    K1 = first(K1),
    K2 = first(K2),
    .groups = "drop"
  )

write.csv(
  trace_summary,
  file.path(out_dir, "summary", "trace_summary_median_l2.csv"),
  row.names = FALSE
)

write.csv(
  budget_lines,
  file.path(out_dir, "summary", "budget_lines.csv"),
  row.names = FALSE
)

# Budget-matched diagnostic table: TS final point vs FD at the corresponding
# dotted vertical budget line. This is only a diagnostic summary; the main plot
# still uses the original iteration axis.
ts_final_by_rep <- raw_results %>%
  filter(arm_type == "TS", eligible == TRUE) %>%
  group_by(rep_id, Config) %>%
  slice_max(order_by = t, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    rep_id,
    Config,
    method = "TS final",
    eval_t = t,
    l2_error
  )

fd_at_budget_by_rep <- raw_results %>%
  filter(arm_type == "TS", eligible == TRUE) %>%
  group_by(rep_id, Config) %>%
  summarise(T_budget = first(T_fd_budget), .groups = "drop") %>%
  left_join(
    raw_results %>%
      filter(arm_type == "FD") %>%
      dplyr::select(rep_id, fd_t = t, fd_l2_error = l2_error),
    by = "rep_id"
  ) %>%
  filter(fd_t == T_budget) %>%
  transmute(
    rep_id,
    Config,
    method = "FD at matched budget",
    eval_t = fd_t,
    l2_error = fd_l2_error
  )

budget_compare_summary <- bind_rows(ts_final_by_rep, fd_at_budget_by_rep) %>%
  group_by(Config, method) %>%
  summarise(
    eval_t = first(eval_t),
    median_l2 = median(l2_error, na.rm = TRUE),
    q25_l2 = quantile(l2_error, 0.25, na.rm = TRUE),
    q75_l2 = quantile(l2_error, 0.75, na.rm = TRUE),
    n_rep = dplyr::n_distinct(rep_id),
    .groups = "drop"
  )

write.csv(
  budget_compare_summary,
  file.path(out_dir, "summary", "budget_matched_TS_final_vs_FD.csv"),
  row.names = FALSE
)

all_configs <- unique(trace_summary$Config)
ts_configs <- sort(all_configs[all_configs != "FD (Baseline)"])

my_colors <- c("FD (Baseline)" = "grey50")
my_linetypes <- c("FD (Baseline)" = "dashed")
my_linewidths <- c("FD (Baseline)" = 1.2)

palette_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                    "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62")

for (i in seq_along(ts_configs)) {
  color_idx <- ((i - 1) %% length(palette_colors)) + 1
  my_colors[ts_configs[i]] <- palette_colors[color_idx]
  my_linetypes[ts_configs[i]] <- "solid"
  my_linewidths[ts_configs[i]] <- 0.8
}

y_max_label <- max(trace_summary$median_l2, na.rm = TRUE)

p_compare_without_title <- ggplot(
  trace_summary,
  aes(x = t, y = median_l2, color = Config, linetype = Config, linewidth = Config)
) +
  geom_line() +
  geom_vline(
    data = budget_lines,
    aes(xintercept = T_budget, color = Config),
    linetype = "dotted",
    linewidth = 1
  ) +
  geom_text(
    data = budget_lines,
    aes(x = T_budget, y = y_max_label, label = paste0("Budget=", T_budget)),
    angle = 90,
    vjust = -0.5,
    hjust = 1,
    size = 3.5,
    show.legend = FALSE
  ) +
  theme_bw() +
  scale_y_log10() +
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(values = my_linetypes) +
  scale_linewidth_manual(values = my_linewidths) +
  labs(
    #title = "Trajectory Comparison: Two-Stage RIGHT vs Full-Data RIGHT",
    #subtitle = "Dotted vertical lines indicate the FD-equivalent sample-access budget for each TS configuration.",
    x = expression(Iteration~Step~(t)),
    y = expression(Median~L[2]~Error),
    color = "Algorithm Configuration",
    linetype = "Algorithm Configuration",
    linewidth = "Algorithm Configuration"
  ) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(3, "lines")
  )

print(p_compare_without_title)



ggsave(
  file.path(out_dir, "figures", "TS_vs_FD_trajectory_plot(no_title).pdf"),
  plot = p_compare_without_title,
  width = 10.5,
  height = 6
)





# ==============================================================================
# Compact table: FD_RIGHT at matched budget vs TS_RIGHT Stage-I/Stage-II endpoint
# ==============================================================================

library(dplyr)

out_dir <- if (exists("out_dir")) {
  out_dir
} else {
  "results/07_schedule_ablation/TS_vs_FD_trajectory"
}

if (!exists("raw_results")) {
  raw_results <- readRDS(file.path(out_dir, "raw", "raw_TS_vs_FD_trajectory.rds"))
}

#dir.create(file.path(out_dir, "summary"), recursive = TRUE, showWarnings = FALSE)

q25 <- function(x) unname(quantile(x, 0.25, na.rm = TRUE))
q75 <- function(x) unname(quantile(x, 0.75, na.rm = TRUE))

# ------------------------------------------------------------------------------
# 1. Extract TS Stage-I endpoint and TS final endpoint, replication by replication
# ------------------------------------------------------------------------------

ts_by_rep <- raw_results %>%
  filter(arm_type == "TS", eligible == TRUE, !is.na(t), !is.na(l2_error)) %>%
  group_by(rep_id, Config) %>%
  summarise(
    q = first(q),
    T1 = first(T1),
    T2 = first(T2),
    B = first(T_fd_budget),
    n1 = first(n1),
    n2 = first(n2),
    b2 = first(b2),
    K1 = first(K1),
    K2 = first(K2),
    
    ts_stage1_l2 = {
      idx <- which(t <= first(T1))
      if (length(idx) == 0L) NA_real_ else l2_error[idx[which.max(t[idx])]]
    },
    
    ts_final_l2 = {
      idx <- which(!is.na(t))
      if (length(idx) == 0L) NA_real_ else l2_error[idx[which.max(t[idx])]]
    },
    
    ts_final_t = {
      idx <- which(!is.na(t))
      if (length(idx) == 0L) NA_real_ else max(t[idx], na.rm = TRUE)
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    ts_stage2_drop_pct = 100 * (1 - ts_final_l2 / ts_stage1_l2)
  )

# ------------------------------------------------------------------------------
# 2. Extract FD error at the matched budget B for each TS configuration
#    If FD trace does not record every integer t, this uses the largest recorded
#    FD t not exceeding B.
# ------------------------------------------------------------------------------

fd_trace <- raw_results %>%
  filter(arm_type == "FD", !is.na(t), !is.na(l2_error)) %>%
  transmute(
    rep_id,
    fd_t = t,
    fd_l2 = l2_error
  )

fd_at_budget_by_rep <- ts_by_rep %>%
  dplyr::select(rep_id, Config, B) %>%
  left_join(fd_trace, by = "rep_id") %>%
  filter(fd_t <= B) %>%
  group_by(rep_id, Config, B) %>%
  slice_max(order_by = fd_t, n = 1, with_ties = FALSE) %>%
  ungroup()

paired_by_rep <- ts_by_rep %>%
  left_join(
    fd_at_budget_by_rep %>%
      dplyr::select(rep_id, Config, B, fd_eval_t = fd_t, fd_l2),
    by = c("rep_id", "Config", "B")
  ) %>%
  mutate(
    fd_over_ts = fd_l2 / ts_final_l2,
    ts_reduction_vs_fd_pct = 100 * (1 - ts_final_l2 / fd_l2)
  )

write.csv(
  paired_by_rep,
  file.path(out_dir, "summary", "compact_table_by_rep.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 3. Numeric summary for checking and manual replacement in LaTeX
# ------------------------------------------------------------------------------

compact_summary <- paired_by_rep %>%
  # Group by core configuration parameters only
  group_by(Config, q, T1, T2) %>%
  summarise(
    # Calculate median L2 errors (removing all q25/q75 calculations as requested)
    fd_median_l2 = median(fd_l2, na.rm = TRUE),
    ts_stage1_median_l2 = median(ts_stage1_l2, na.rm = TRUE),
    ts_final_median_l2 = median(ts_final_l2, na.rm = TRUE),
    
    # Calculate the median percentage drop for Stage II
    median_stage2_drop_pct = median(ts_stage2_drop_pct, na.rm = TRUE),
    
    # Optional relative metrics for sanity check
    median_fd_over_ts = median(fd_over_ts, na.rm = TRUE),
    median_ts_reduction_vs_fd_pct = median(ts_reduction_vs_fd_pct, na.rm = TRUE),
    
    n_rep = n_distinct(rep_id),
    .groups = "drop"
  ) %>%
  # Sort by hyper-parameters to ensure chronological tabular output
  arrange(T1, T2, q)

write.csv(
  compact_summary,
  file.path(out_dir, "summary", "compact_table_summary_numeric.csv"),
  row.names = FALSE
)

print(compact_summary)


# ------------------------------------------------------------------------------
# 4. Generate a compact LaTeX tabular fragment (Column-based layout)
# ------------------------------------------------------------------------------

digits <- 3

# Helper function: Format numbers with fixed digits, return "--" for NAs
fmt_num <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "--",
    formatC(x, format = "f", digits = digits)
  )
}

# Helper function: Wrap number in LaTeX math mode, apply bold conditionally
fmt_error <- function(m, bold = FALSE, digits = 3) {
  m_s <- fmt_num(m, digits)
  if (bold) {
    sprintf("$\\mathbf{%s}$", m_s)
  } else {
    sprintf("$%s$", m_s)
  }
}

# Helper function: Format the percentage drop/increase with arrows
fmt_drop <- function(x) {
  if (is.na(x)) return("")
  if (x >= 0) {
    sprintf("$\\downarrow %.1f\\%%$", x)
  } else {
    sprintf("$\\uparrow %.1f\\%%$", abs(x))
  }
}

# Generate row strings for the LaTeX table
latex_rows <- compact_summary %>%
  rowwise() %>%
  transmute(
    # Column 1: Method Configuration
    Method = sprintf(
      "$\\mathrm{TS}_{\\mathrm{RIGHT}}$ ($q=%.2f$, $(T_1,T_2)=(%d,%d)$)",
      q, T1, T2
    ),
    
    # Column 2: FD final error (Baseline placed in a column instead of a row)
    FD_Final = fmt_error(fd_median_l2, bold = FALSE, digits = digits),
    
    # Column 3: TS Stage I error
    TS_Stage1 = fmt_error(ts_stage1_median_l2, bold = FALSE, digits = digits),
    
    # Column 4: TS Final error with stacked percentage drop
    TS_Final = sprintf(
      "\\shortstack[c]{%s\\\\{\\scriptsize %s from Stage I}}",
      fmt_error(ts_final_median_l2, bold = TRUE, digits = digits),
      fmt_drop(median_stage2_drop_pct)
    )
  ) %>%
  ungroup() %>%
  # Concatenate columns with ' & ' and append '\\' for LaTeX newline
  mutate(
    latex = paste(Method, FD_Final, TS_Stage1, TS_Final, sep = " & "),
    latex = paste0(latex, " \\\\")
  )

# Construct the complete LaTeX table environment
latex_fragment <- c(
  "% Required packages: booktabs, array",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{4pt}",       # Slightly increased padding for readability
  "\\renewcommand{\\arraystretch}{1.2}", # Slightly increased height to accommodate \shortstack
  "\\begin{tabular}{lccc}",              # 4 columns: Left, Center, Center, Center
  "\\toprule",
  # Updated Header: FD is now explicitly labeled as a column
  "Method & $\\mathrm{FD}_{\\mathrm{RIGHT}}$ (Matched Budget) & Stage-I endpoint & Final / Stage-II endpoint \\\\",
  "\\midrule",
  latex_rows$latex,
  "\\bottomrule",
  "\\end{tabular}"
)

writeLines(
  latex_fragment,
  file.path(out_dir, "summary", "compact_budget_table_fragment.tex")
)

# Print output to console for immediate preview
cat(paste(latex_fragment, collapse = "\n"))
