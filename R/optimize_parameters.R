
#' optimize_parameters
#'
#' Optimize parameters for the biphasic exponential decay model using non-linear least squares.
#' 
#' @param data prepared data with required time and F_t columns
#' @param initial_params a vector of numeric values; for biphasic exponential decay model input c(F, S, k_f, k_s)
#'
#' @return optimized parameters in c(F, S, k_f, k_s)
#' @export

optimize_parameters <- function(data, initial_params) {
  model <- function(t, F, S, k_f, k_s) {
    F * exp(-k_f * t) + S * exp(-k_s * t)
  }
  # Define the objective function internally
  objective <- function(params) {
    F <- params[1]
    S <- params[2]
    k_f <- params[3]
    k_s <- params[4]
    predicted <- model(data$time, F, S, k_f, k_s)
    sum((data$F_t - predicted)^2)  # Sum of squared differences
  }

  # Call optim to optimize parameters
  tryCatch({
    opt_result <- optim(initial_params, objective)
    return(opt_result$par)
  }, error = function(e) {
    # Print error message and exit function
    message("Error in optim: ", e$message)
    return(NULL)
  })
}


