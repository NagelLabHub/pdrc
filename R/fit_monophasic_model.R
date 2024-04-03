# 8. monophasic decay modeling

#' fit_monophasic_model
#'
#' @import minpack.lm
#' @import ggthemes
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise n
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_line coord_cartesian scale_y_continuous annotate xlab ylab
#' @importFrom stats coef optim predict residuals sd time uniroot
#'
#' @param data prepared data for monophasic model
#' @param initial_params a vector of two values c(F, k)
#'
#' @return list of results
#' @export
#'
#' @examples
#' out_data <- data.frame(time = c(0, 15, 30, 60, 120), F_t = c(48,40,19,15,12))
#' mono_list <- fit_monophasic_model(out_data, list(F = 48, k = 0.04))
fit_monophasic_model <- function(data, initial_params) {
  # Define the monophasic model function
  model_monophasic <- function(time, F, k) {
    F * exp(-k * time)
  }

  # Fit the monophasic model using provided initial parameters
  fit_monophasic <- minpack.lm::nlsLM(F_t ~ model_monophasic(time, F, k),
                        data = data,
                        start = initial_params)

  sample_name <- unique(data$Sample_Name)

  # Extract fitted parameters for monophasic model
  F_fit_monophasic <- coef(fit_monophasic)[1]
  k_fit_monophasic <- coef(fit_monophasic)[2]

  # Calculate half-life for monophasic model
  half_life_overall_monophasic <- log(2) / k_fit_monophasic

  # Calculate statistics for monophasic model
  residuals_monophasic <- residuals(fit_monophasic)
  n <- length(data$F_t)  # Number of observations
  p_monophasic <- length(coef(fit_monophasic))  # Number of parameters
  RSE_monophasic <- sqrt(sum(residuals_monophasic^2) / (n - p_monophasic))
  SS_total_monophasic <- sum((data$F_t - mean(data$F_t))^2)  # Total sum of squares
  SS_residual_monophasic <- sum(residuals_monophasic^2)  # Residual sum of squares
  R_squared_monophasic <- 1 - (SS_residual_monophasic / SS_total_monophasic)

  new_data <- data.frame(time = seq(min(data$time), max(data$time), length.out = 1000))
  new_data$F_t_predicted_monophasic <- predict(fit_monophasic, newdata = new_data)

  summary_data <- data %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(avg_F_t = mean(F_t),
              sem_F_t = sd(F_t) / sqrt(dplyr::n()))

  # Create ggplot object for the monophasic model
  plot_object_monophasic <- ggplot() +
    geom_point(data = summary_data, aes(x = time, avg_F_t), color = "blue", size=2) +
    geom_errorbar(data = summary_data, aes(x = time, ymin = avg_F_t - sem_F_t, ymax = avg_F_t + sem_F_t), color = "blue", size = 0.8, width=2) +
    ggthemes::theme_few() +
    geom_line(data = new_data, aes(x = time, y = F_t_predicted_monophasic), linetype = "dashed") + coord_cartesian(clip = "off") +
    scale_y_continuous(limits = c(0, 60)) +
    annotate("text", x = Inf, y = Inf, label = sprintf("Sample: %s\nR_sqaured = %.2f\nt_half overall = %.2f", sample_name, R_squared_monophasic, half_life_overall_monophasic), hjust = 1.1, vjust = 2, size = 4, colour = "black")  +
    xlab("Repair time (min)") +
    ylab("% DNA in tail \n (background corrected)")

  # Return the fitted parameters, half-life, RSE, R-squared, and plot object for the monophasic model
  return(list(
    fit = fit_monophasic,
    params = c(F_fit_monophasic, k_fit_monophasic),
    half_life = half_life_overall_monophasic,
    RSE = RSE_monophasic,
    R_squared = R_squared_monophasic,
    plot_object = plot_object_monophasic,
    sample_name = sample_name
  ))
}

loop_data_mono <- function(dataset_name) {

  results_list <- list()

  for (sample_name in unique(dataset_name$Sample)) {
    data <- prepare_data(dataset_name, sample_name)

    initial_params_list <- list(
      c(F = data$F_t[1], k = 0.04)
    )

    fitted <- FALSE
    for (initial_params in initial_params_list) {
      tryCatch({
        fit_results <- fit_monophasic_model(data, initial_params)

        results_list[[sample_name]] <- fit_results

        fitted <- TRUE
        break  # Exit the loop if successful fit
      }, error = function(e) {
        # Print error message and try next set of initial_params
        cat(paste("Error processing sample:", sample_name, "- Trying next set of initial parameters.\n"))
      })

      if (fitted) break  # Exit the loop if fitted successfully
    }

    if (!fitted) {
      # If none of the initial_params combinations worked, add an empty row
      cat(paste("Error processing sample:", sample_name, "- None of the initial parameters combinations worked.\n"))
    }
  }

  # Return the results_list
  return(results_list)
}

convert_results_to_dataframe_mono <- function(results_list) {
  # Initialize an empty dataframe to store the results
  list_dat <- data.frame(Sample = character(),
                         Coef_F = numeric(),
                         Coef_k = numeric(),
                         R_squared = numeric(),
                         Half_life_overall = numeric(),
                         stringsAsFactors = FALSE)

  # Iterate over each element in results_list
  for (sample_name in names(results_list)) {
    sample_results <- results_list[[sample_name]]  # Extract results for the current sample

    # Create a row for the current sample and append it to the list_dat dataframe
    list_dat <- rbind(list_dat, data.frame(Sample = sample_name,
                                           Coef_F = sample_results$params[1],
                                           Coef_k = sample_results$params[2],
                                           R_squared = sample_results$R_squared,
                                           Half_life_overall = sample_results$half_life,
                                           stringsAsFactors = FALSE))
  }

  # Return the prebatch dataframe
  return(list_dat)
}
