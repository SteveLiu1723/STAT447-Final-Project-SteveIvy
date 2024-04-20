compute_parameter_mcse <- function(mcmc_list) {
  calc_mcse <- function(mcmc_obj) {
    last_samples <- window(mcmc_obj, start = end(mcmc_obj) - 9999)
    mcse_vals <- summary(last_samples)$statistics[, "Time-series SE"]
    return(mcse_vals)
  }
  
  mcse_matrix <- t(sapply(mcmc_list, calc_mcse))
  avg_mcse_per_param <- colMeans(mcse_matrix, na.rm = TRUE)
  return(avg_mcse_per_param)
}