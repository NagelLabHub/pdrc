#' setup_parallel
#'
#' Setup parallelization for faster Bayesian model fitting and adaptive spike modeling.
#'
#' @param workers Number of CPU cores to use. Default = all available cores.
#'
#' @return None (sets options globally).
#' @export
#' @importFrom future plan multisession availableCores
#' @importFrom rstan rstan_options
setup_parallel <- function(workers = future::availableCores()) {
  library(future)
  library(rstan)

  plan(multisession, workers = workers)
  options(mc.cores = workers)
  options(brms.backend = "rstan")
  rstan_options(auto_write = TRUE, javascript = FALSE)
  Sys.setenv(STAN_NUM_THREADS = workers)
  Sys.setenv(RCPP_PARALLEL_BACKEND = "tbb")

  message(sprintf("Parallelization successfully set up with %d workers.", workers))
}
