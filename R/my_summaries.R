#' @import minpack.lm
#' @import tidyverse

# 1. biphasic exponential decay model

model <- function(t, F, S, k_f, k_s) {
  F * exp(-k_f * t) + S * exp(-k_s * t)
}

# 2. data prep - separate for replicates input vs averaged input

prepare_data <- function(dataset, sample_name) {
  time <- c(0, 15, 30, 60, 120)

  columns_to_use <-c("c_0", "c_15", "c_30", "c_60", "c_120")

  F_t <- dataset[dataset$Sample == sample_name, columns_to_use]

  data <- data.frame(
    time = rep(time, nrow(F_t)),
    F_t = as.vector(t(F_t)),
    Sample_Name = rep(sample_name, length(time) * nrow(F_t))
  )

  return(data)
}

# 3. fine tuning for parameter estimation

## Define the objective function for optimization
objective <- function(params, data, model_function) {
  F <- params[1]
  S <- params[2]
  k_f <- params[3]
  k_s <- params[4]
  predicted <- model_function(data$time, F, S, k_f, k_s)
  sum((data$F_t - predicted)^2)  # Sum of squared differences
}

## Function to perform optimization and return optimized parameters
optimize_parameters <- function(data, model_function, initial_params) {
  opt_result <- optim(initial_params, objective, data = data, model_function = model_function)
  return(opt_result$par)
}

# 4. half-life estimation for overall, fast, slow phases

calculate_half_life <- function(F, S, k_f, k_s) {
  half_life_eq <- function(t_half, F, S, k_f, k_s) {
    F * exp(-k_f * t_half) + S * exp(-k_s * t_half) - (F + S) / 2 }

  half_life_eq_fast <- function(t_half, k) {
    exp(-k * t_half) - 0.5 }

  half_life_eq_slow <- function(t_half, k) {
    exp(-k * t_half) - 0.5 }

  initial_interval <- c(0, 360)
  slow_phase_result <- 360

  # Use tryCatch to handle potential errors
  tryCatch({

    root <- uniroot(half_life_eq, interval = initial_interval, F = F, S = S, k_f = k_f, k_s = k_s)
    root_fast <- uniroot(half_life_eq_fast, interval = initial_interval, k = k_f)
    root_slow <- uniroot(half_life_eq_slow, interval = initial_interval, k = k_s)

    # Check if slow phase root was successful
    if (!is.null(root_slow$root)) {
      slow_phase_result <- root_slow$root
    }
  }, error = function(e) {
    # Do nothing or handle the error as needed
    # Omitting error output in this case
  })

  return(list(overall = ifelse(is.null(root$root), NA, root$root),
              fast_phase = ifelse(is.null(root_fast$root), NA, root_fast$root),
              slow_phase = slow_phase_result))
}

# 5. main model fitting process

fit_model_and_stats <- function(data, model_function, initial_params, calculate_half_life) {
  # Fit the model using provided initial parameters
  fit <- nlsLM(F_t ~ model_function(time, F, S, k_f, k_s),
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
    group_by(time) %>%
    summarise(avg_F_t = mean(F_t),
              sem_F_t = sd(F_t) / sqrt(n()))

  # Create a ggplot object for the data and the model prediction
  plot_object <- ggplot() +
    geom_point(data = summary_data, aes(x = time, avg_F_t), color = "blue", size=2) +
    geom_errorbar(data = summary_data, aes(x=time, ymin=avg_F_t-sem_F_t, ymax=avg_F_t+sem_F_t), color = "blue", size = 0.8, width=2) +
    theme_few() +
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

# 6. loop through all samples and apply the model fitting in step 5.

loop_data <- function(dataset_name) {

  results_list <- list()

  for (sample_name in unique(dataset_name$Sample)) {
    data <- prepare_data(dataset_name, sample_name)

    initial_params_list <- list(
      c(F = data$F_t[1], S = 0, k_f = 0.04, k_s = 0.02),
      optimize_parameters(data, model, c(F = data$F_t[1], S = 0, k_f = 0.04, k_s = 0.02)),
      c(F = data$F_t[1], S = 0, k_f = 0.4, k_s = 0.04),
      c(F = 0.5 * data$F_t[1], S = 0.5 * data$F_t[1], k_f = 0.4, k_s = 0.04),
      c(F = data$F_t[1], S = 0, k_f = 0.1, k_s = 0.01),
      c(F = 0.5 * data$F_t[1], S = 0.5 * data$F_t[1], k_f = 0.1, k_s = 0.01),
      c(F = data$F_t[1], S = 0, k_f = 0.04, k_s = 0)
    )

    fitted <- FALSE
    for (initial_params in initial_params_list) {
      tryCatch({
        fit_results <- fit_model_and_stats(data, model, initial_params, calculate_half_life)

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

# 7. summarize params into dataframe from returned list in step 6.

convert_results_to_dataframe <- function(results_list) {
  # Initialize an empty dataframe to store the results
  list_dat <- data.frame(Sample = character(),
                         Coef_F = numeric(),
                         Coef_S = numeric(),
                         Coef_k_f = numeric(),
                         Coef_k_s = numeric(),
                         R_squared = numeric(),
                         Half_life_overall = numeric(),
                         Half_life_fast = numeric(),
                         Half_life_slow = numeric(),
                         stringsAsFactors = FALSE)

  # Iterate over each element in results_list
  for (sample_name in names(results_list)) {
    sample_results <- results_list[[sample_name]]  # Extract results for the current sample

    # Create a row for the current sample and append it to the list_dat dataframe
    list_dat <- rbind(list_dat, data.frame(Sample = sample_name,
                                           Coef_F = sample_results$params[1],
                                           Coef_S = sample_results$params[2],
                                           Coef_k_f = sample_results$params[3],
                                           Coef_k_s = sample_results$params[4],
                                           R_squared = sample_results$R_squared,
                                           Half_life_overall = sample_results$half_life$overall,
                                           Half_life_fast = sample_results$half_life$fast_phase,
                                           Half_life_slow = sample_results$half_life$slow_phase,
                                           stringsAsFactors = FALSE))
  }

  # Return the prebatch dataframe
  return(list_dat)
}

# 8. monophasic decay modeling
fit_monophasic_model <- function(data, initial_params) {
  # Define the monophasic model function
  model_monophasic <- function(time, F, k) {
    F * exp(-k * time)
  }

  # Fit the monophasic model using provided initial parameters
  fit_monophasic <- nls(F_t ~ model_monophasic(time, F, k),
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
    group_by(time) %>%
    summarise(avg_F_t = mean(F_t),
              sem_F_t = sd(F_t) / sqrt(n()))

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
