# 5. main model fitting process

#' fit_model_and_stats
#'
#' @import minpack.lm
#' @import ggthemes
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise n
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_line coord_cartesian scale_y_continuous annotate xlab ylab
#' @importFrom stats coef optim predict residuals sd time uniroot
#'
#' @param data prepared data
#' @param model_function model formula
#' @param initial_params initial values for fitting
#'
#' @return a list of model fit par, estimates and plot
#' @export
#'
#' @examples
#' out_data <- data.frame(time = c(0, 15, 30, 60, 120), F_t = c(40,26,19,15,12))
#' params <- c(40,0,0.1,0.01)
#' result_list = fit_model_and_stats(out_data, model, params)
#'
fit_model_and_stats <- function(data, model_function, initial_params) {
  # Fit the model using provided initial parameters
  fit <- minpack.lm::nlsLM(F_t ~ model_function(time, F, S, k_f, k_s),
               data = data, control = list(maxiter = 200),
               lower = c(F = 0, S = 0, k_f = 0, k_s = 0),
               start = list(F = initial_params[1],
                            S = initial_params[2],
                            k_f = initial_params[3],
                            k_s = initial_params[4]))

  sample_name <- unique(data$Sample_Name)

  # Extract fitted parameters
  F_fit <- coef(fit)[1]
  S_fit <- coef(fit)[2]
  k_f_fit <- coef(fit)[3]
  k_s_fit <- coef(fit)[4]

  # Check if k_s > k_f and switch if necessary
  if (k_s_fit > k_f_fit) {
    temp_k <- k_s_fit
    k_s_fit <- k_f_fit
    k_f_fit <- temp_k

    temp_S <- S_fit
    S_fit <- F_fit
    F_fit <- temp_S
  }

  half_life <- calculate_half_life(F_fit, S_fit, k_f_fit, k_s_fit)

  half_life_overall <- half_life$overall
  half_life_fast <- half_life$fast_phase
  half_life_slow <- half_life$slow_phase

  # Calculate statistics
  residuals <- residuals(fit)
  n <- length(data$F_t)  # Number of observations
  p <- length(coef(fit))  # Number of parameters
  RSE <- sqrt(sum(residuals^2) / (n - p))
  SS_total <- sum((data$F_t - mean(data$F_t))^2)  # Total sum of squares
  SS_residual <- sum(residuals^2)  # Residual sum of squares
  R_squared <- 1 - (SS_residual / SS_total)

  new_data <- data.frame(time = seq(min(data$time), max(data$time), length.out = 1000))
  new_data$F_t_predicted <- predict(fit, newdata = new_data)

  summary_data <- data %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(avg_F_t = mean(F_t),
              sem_F_t = sd(F_t) / sqrt(n()))

  # Create a ggplot object for the data and the model prediction
  plot_object <- ggplot() +
    geom_point(data = summary_data, aes(x = time, avg_F_t), color = "blue", size=2) +
    geom_errorbar(data = summary_data, aes(x=time, ymin=avg_F_t-sem_F_t, ymax=avg_F_t+sem_F_t), color = "blue", size = 0.8, width=2) +
    ggthemes::theme_few() +
    geom_line(data = new_data, aes(x = time, y = F_t_predicted), linetype = "dashed") + coord_cartesian(clip = "off") +
    scale_y_continuous(limits = c(0, 60)) +
    annotate("text", x = Inf, y = Inf, label = sprintf("Sample: %s\nR_sqaured = %.2f\nt_half overall = %.2f\nt_half fast = %.2f", sample_name, R_squared, half_life_overall, half_life_fast), hjust = 1.1, vjust = 2, size = 4, colour = "black")  +
    xlab("Repair time (min)") +
    ylab("% DNA in tail \n (background corrected)")

  # Return the fitted parameters, half-life, RSE, R-squared, and plot object
  return(list(fit = fit,
              params = c(F_fit, S_fit, k_f_fit, k_s_fit),
              half_life = half_life,
              RSE = RSE,
              R_squared = R_squared,
              plot_object = plot_object,
              sample_name = sample_name))
}
