script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  }
  getwd()
}
wd0 <- getwd()
setwd(script_dir())
on.exit(setwd(wd0), add = TRUE)

# Convenience driver. Running all diagnostics can be slow.
source("31_support_localization_FD_linear_parallel.R")
source("32_support_localization_FD_logistic_parallel.R")
source("33_gradient_complexity_pointwise_linear_parallel.R")
source("34_gradient_complexity_pointwise_logistic_parallel.R")
source("35_gradient_complexity_uniform_linear_parallel.R")
source("36_gradient_complexity_uniform_logistic_parallel.R")
