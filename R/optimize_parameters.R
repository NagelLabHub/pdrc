# 3. fine tuning for parameter estimation

#' optimize_parameters
#'
#' @param data prepared data with required time and F_t columns
#' @param model_function model formula to be optimized
#' @param initial_params a vector of numeric values; for biphasic exponential decay model input c(F, S, k_f, k_s)
#'
#' @return optimized parameters in c(F, S, k_f, k_s)
#' @export
#'
#' @examples
#' out_data <- data.frame(time = c(0, 15, 30, 60, 120), F_t = c(40,26,19,15,12))
#' params <- optimize_parameters(out_data, model, c(40,0,0.1,0.01))
optimize_parameters <- function(data, model_function, initial_params) {
  # Define the objective function internally
  objective <- function(params) {
    F <- params[1]
    S <- params[2]
    k_f <- params[3]
    k_s <- params[4]
    predicted <- model_function(data$time, F, S, k_f, k_s)
    sum((data$F_t - predicted)^2)  # Sum of squared differences
  }

  # Call optim to optimize parameters
  opt_result <- optim(initial_params, objective)
  return(opt_result$par)
}

