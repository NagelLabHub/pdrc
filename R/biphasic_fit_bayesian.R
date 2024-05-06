# 4. main model fitting process - Bayesian inference

#' biphasic_fit_bayesian
#'
#' @import brms
#' @import ggthemes
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise n
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_line coord_cartesian scale_y_continuous annotate xlab ylab
#' @importFrom stats coef optim predict residuals sd time uniroot
#'
#' @param data dataset prepared by prepare_data
#'
#' @return a large list of result lists ordered by sample name
#' @export
#'
#' @examples
#' single_data <- data.frame(time = c(0, 15, 30, 60, 120), F_t = c(40,26,19,15,12))
#' result_list <- biphasic_fit_bayesian(single_data)
#'
biphasic_fit_bayesian <- function(data) {
  prior <- c(
  set_prior("normal(30, 10)", nlpar = "F", lb = 0),  # F is greater than or equal to 0
  set_prior("normal(10, 10)", nlpar = "S", lb = 0),  # S is greater than or equal to 0
  set_prior("normal(0.04, 0.04)", nlpar = "kf", lb = 0),  # k_f is greater than 0
  set_prior("normal(0.01, 0.01)", nlpar = "ks", lb = 0) )  # k_s is greater than or equal to 0

  formula <- bf(F_t ~ F * exp(-kf * time) + S * exp(-ks * time),
              F ~ 1,  # non-linear parameter
              S ~ 1,  # non-linear parameter
              kf ~ 1,  # non-linear parameter
              ks ~ 1,  # non-linear parameter
              nl = TRUE)

  b_fit <- brm(formula,
                 data = data,
                 prior = prior,
                 control = list(adapt_delta = 0.95, max_treedepth = 15),
                 chains = 4,
                 iter = 2000,
                 warmup = 1000)

# Extract fitted parameters
F_fit <- fixef(b_fit)[1,1]
S_fit <- fixef(b_fit)[2,1]
k_f_fit <- fixef(b_fit)[3,1]
k_s_fit <- fixef(b_fit)[4,1]

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

loo <- loo(b_fit)
looic <- loo$estimates["looic","Estimate"]

new_data <- data.frame(time = seq(min(data$time), max(data$time), length.out = 1000))
new_data$F_t_predicted <- predict(b_fit, newdata = new_data)[,1]

summary_data <- data %>%
  group_by(time) %>%
  summarise(avg_F_t = mean(F_t),
            sem_F_t = sd(F_t) / sqrt(n()))

sample_name <- unique(data$Sample_Name)

plot_object <- ggplot() +
  geom_point(data = summary_data, aes(x = time, avg_F_t), color = "blue", size=2) +
  geom_errorbar(data = summary_data, aes(x=time, ymin=avg_F_t-sem_F_t, ymax=avg_F_t+sem_F_t), color = "blue", linewidth = 0.8, width=2) +
  ggthemes::theme_few() +
  geom_line(data = new_data, aes(x = time, y = F_t_predicted), linetype = "dashed") + coord_cartesian(clip = "off") +
  scale_y_continuous(limits = c(0, 60)) +
  annotate("text", x = Inf, y = Inf, label = sprintf("Sample: %s\nLOOIC = %.2f\nt_half overall = %.2f\nt_half fast = %.2f\nt_half slow = %.2f", sample_name, looic, half_life_overall, half_life_fast, half_life_slow), hjust = 1.1, vjust = 2, size = 4, colour = "black")  +
  xlab("Repair time (min)") +
  ylab("% DNA in tail \n (background corrected)")

# Return the fitted parameters, half-life, RSE, R-squared, and plot object
return(list(b_fit = b_fit,
            params = c(F_fit, S_fit, k_f_fit, k_s_fit),
            half_life = half_life,
            loo = loo,
            plot_object = plot_object,
            sample_name = sample_name))
}
