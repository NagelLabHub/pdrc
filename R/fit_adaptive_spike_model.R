#' fit_adaptive_spike_model
#'
#' Fit multiple decay models (fast decay, slow decay, accumulation) and select the best model based on LOOIC, spike detection, biological penalties, and model complexity.
#'
#' @param data A data frame with 'time' and 'F_t' columns.
#' @param bio_penalty Numeric; penalty applied for biologically implausible models (default = 3).
#' @param complexity_penalty Numeric; penalty per additional model term (default = 2).
#'
#' @return A list containing the best model, model scores, and flags for spike detection and short follow-up.
#' @export
#' 
#' @importFrom brms brm bf set_prior loo
#' @importFrom dplyr arrange pull mutate filter
#' @importFrom furrr future_map furrr_options
#' @importFrom purrr compact

fit_adaptive_spike_model <- function(data, bio_penalty = 4, complexity_penalty = 2) {
  
  sample_name <- unique(data$Sample_Name)
  
  # Detect spike
  time_sorted <- data %>% arrange(time)
  t0_value <- mean(time_sorted$F_t[time_sorted$time == min(time_sorted$time)])
  next_time <- sort(unique(time_sorted$time))[2]
  t_next_value <- mean(time_sorted$F_t[time_sorted$time == next_time])
  spike_detected <- (t_next_value - t0_value) / t0_value > 0.05
  
  # Detect follow-up duration
  max_time <- max(data$time)
  short_followup <- max_time <= 60
  
  message(sprintf("Spike detected: %s", spike_detected))
  message(sprintf("Short follow-up (<=60 min): %s", short_followup))
  
  # Define models
  formula_fast <- bf(F_t ~ F * exp(-k2 * time), F + k2 ~ 1, nl = TRUE)
  formula_fast_slow <- bf(F_t ~ F * exp(-k2 * time) + S * exp(-k3 * time), F + S + k2 + k3 ~ 1, nl = TRUE)
  formula_fast_accum <- bf(F_t ~ F * exp(-k2 * time) + A * (1 - exp(-k1 * time)), F + A + k1 + k2 ~ 1, nl = TRUE)
  formula_fast_slow_accum <- bf(F_t ~ F * exp(-k2 * time) + S * exp(-k3 * time) + A * (1 - exp(-k1 * time)),
                                F + S + A + k1 + k2 + k3 ~ 1, nl = TRUE)
  
  # Priors
  priors_fast <- c(
    set_prior("normal(40, 10)", nlpar = "F", lb = 0),
    set_prior("normal(0.2, 0.1)", nlpar = "k2", lb = 0)
  )
  priors_fast_slow <- c(
    set_prior("normal(40, 10)", nlpar = "F", lb = 0),
    set_prior("normal(10, 10)", nlpar = "S", lb = 0),
    set_prior("normal(0.02, 0.01)", nlpar = "k3", lb = 0),
    set_prior("normal(0.2, 0.1)", nlpar = "k2", lb = 0)
  )
  priors_fast_accum <- c(
    set_prior("normal(40, 10)", nlpar = "F", lb = 0),
    set_prior("normal(0, 10)", nlpar = "A", lb = 0),
    set_prior("normal(0.4, 0.2)", nlpar = "k1", lb = 0),
    set_prior("normal(0.2, 0.1)", nlpar = "k2", lb = 0)
  )
  priors_fast_slow_accum <- c(
    set_prior("normal(40, 10)", nlpar = "F", lb = 0),
    set_prior("normal(10, 10)", nlpar = "S", lb = 0),
    set_prior("normal(0, 10)", nlpar = "A", lb = 0),
    set_prior("normal(0.4, 0.2)", nlpar = "k1", lb = 0),
    set_prior("normal(0.02, 0.01)", nlpar = "k3", lb = 0),
    set_prior("normal(0.2, 0.1)", nlpar = "k2", lb = 0)
  )
  
  model_list <- list(
    fast = list(formula = formula_fast, priors = priors_fast),
    fast_slow = list(formula = formula_fast_slow, priors = priors_fast_slow),
    fast_accum = list(formula = formula_fast_accum, priors = priors_fast_accum),
    fast_slow_accum = list(formula = formula_fast_slow_accum, priors = priors_fast_slow_accum)
  )
  
  # Fitting function
  fit_model <- function(formula, priors) {
    suppressMessages(
      suppressWarnings(
        brm(
          formula = formula, data = data, prior = priors,
          chains = 4, iter = 4000, warmup = 2000, refresh = 0,
          control = list(adapt_delta = 0.99, max_treedepth = 15),
          silent = TRUE
        )
      )
    )
  }
  
  fits <- furrr::future_map(model_list, ~ fit_model(.x$formula, .x$priors), .options = furrr_options(seed = TRUE))
  names(fits) <- names(model_list)
  
  # Calculate LOOICs
  looics <- furrr::future_map(fits, ~ loo(.x, moment_match = TRUE, reloo = TRUE))
  looic_values <- sapply(looics, function(x) x$estimates["looic", "Estimate"])
  model_terms <- c(fast = 1, fast_slow = 2, fast_accum = 2, fast_slow_accum = 3)
  
  # Build score table
  model_scores <- tibble(
    Model = names(looic_values),
    LOOIC = looic_values,
    BioPenalty = 0,
    ComplexityPenalty = 0,
    FinalScore = looic_values,
    Selected = "No"
  )
  
  # Apply soft penalties
  for (i in seq_len(nrow(model_scores))) {
    model <- model_scores$Model[i]

    # Apply biological penalties
    if (!spike_detected && grepl("accum", model)) {
      model_scores$BioPenalty[i] <- model_scores$BioPenalty[i] + bio_penalty
    }
    if (short_followup && grepl("slow", model)) {
      model_scores$BioPenalty[i] <- model_scores$BioPenalty[i] + bio_penalty
    }
    if (spike_detected && grepl("^fast$", model)) {
      model_scores$BioPenalty[i] <- model_scores$BioPenalty[i] + bio_penalty
    }

    # Complexity penalty
    if (model != "fast") {
      model_scores$ComplexityPenalty[i] <- complexity_penalty * (model_terms[model] - 1)
    }

    # Final Score
    model_scores$FinalScore[i] <- model_scores$LOOIC[i] +
    model_scores$BioPenalty[i] +
    model_scores$ComplexityPenalty[i]
  }
  
  # Select model
  best_model_name <- model_scores %>%
    arrange(FinalScore) %>%
    slice(1) %>%
    pull(Model)
  
  best_model <- fits[[best_model_name]]
  
  model_scores$Selected <- ifelse(model_scores$Model == best_model_name, "Yes", "No")
  
  message(sprintf("Selected model: %s (LOOIC = %.2f)", best_model_name, looic_values[best_model_name]))
  
  return(list(
    best_model = best_model,
    best_model_name = best_model_name,
    model_scores = model_scores,
    spike_detected = spike_detected,
    short_followup = short_followup,
    data = data 
  ))
}
