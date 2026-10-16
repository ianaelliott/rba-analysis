#------------------------------------------------------------------------------#
## EO INDEX EXPLORATION v1 ##
#------------------------------------------------------------------------------#

library(stats)
library(ggpubr)
library(tidyverse)

# Set the back-end
#rm(list=ls()) # Clear environment
options(scipen = 999) # Disable scientific notation
seed1 <- 6012021
seed2 <- 3222023
seed3 <- 9061991

# Set the parameters
base_rates <- seq(0.05, 0.85, by = 0.05)
err_rate   <- 0.15   # held constant
n <- 12000

# Apply the base rates
eo <- sapply(base_rates, function(b) {
  score <- rbeta(n, b * 10, (1 - b) * 10)
  obs   <- rbinom(n, 1, score)
  flip  <- runif(n) < err_rate
  obs_e <- ifelse(flip, 1 - obs, obs)
  sum(score) / sum(obs_e)
})

# Bootstrap the EO function
n_boot <- 2000
set.seed(seed1)
boot_curve <- lapply(base_rates, function(b) {
  score <- rbeta(n, b * 10, (1 - b) * 10)
  obs   <- rbinom(n, 1, score)
  replicate(n_boot, {
    idx  <- sample(n, replace = TRUE)
    fl_b <- runif(length(idx)) < err_rate
    sum(score[idx]) / sum(ifelse(fl_b, 1 - obs[idx], obs[idx]))
  })
})

# Collate the outcomes
boot_summary <- do.call(rbind, lapply(seq_along(base_rates), function(i) {
  x <- boot_curve[[i]]
  data.frame(base_rate = base_rates[i],
             eo  = mean(x),
             lwr = quantile(x, .025),
             upr = quantile(x, .975))
}))
print(boot_summary, row.names = FALSE)

# Plot the EO effect
eo_effect_plot <- ggplot(boot_summary, 
       aes(base_rate, eo)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "steelblue") +
  labs(x = "Base rate", y = "Expected/Observed Index", 
       title = "Bootstrapped EO Indices at various base rates (Error rate = 0.15)") +
  theme_bw()
print(eo_effect_plot)

# Brier version
set.seed(seed3)
boot_brier <- lapply(base_rates, function(b) {
  score <- rbeta(n, b * 10, (1 - b) * 10)
  obs   <- rbinom(n, 1, score)
  cal <- err <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    idx   <- sample(n, replace = TRUE)
    fl_b  <- runif(length(idx)) < err_rate
    obs_b <- ifelse(fl_b, 1 - obs[idx], obs[idx])
    cal[i] <- mean((score[idx] - obs[idx])^2)
    err[i] <- mean((score[idx] - obs_b)^2)
  }
  data.frame(base_rate  = b,
             brier_cal = mean(cal), cal_lwr = quantile(cal, .025), cal_upr = quantile(cal, .975),
             brier_err = mean(err), err_lwr = quantile(err, .025), err_upr = quantile(err, .975))
})
boot_brier <- do.call(rbind, boot_brier)
print(boot_brier, row.names = FALSE)

# Plot the Brier effect
brier_effect_plot <- ggplot(boot_brier, aes(base_rate)) +
  geom_ribbon(aes(ymin = cal_lwr, ymax = cal_upr), alpha = 0.15, fill = "darkred") +
  geom_ribbon(aes(ymin = err_lwr, ymax = err_upr), alpha = 0.20, fill = "steelblue") +
  geom_line(aes(y = brier_cal), colour = "darkred") +
  geom_line(aes(y = brier_err), colour = "steelblue") +
  geom_point(aes(y = brier_cal), colour = "darkred") +
  geom_point(aes(y = brier_err), colour = "steelblue") +
  labs(x = "Base rate", y = "Brier score") +
  theme_bw()
print(brier_effect_plot)

brier_effect_plot2 <- ggplot(boot_brier, aes(base_rate)) +
  geom_ribbon(aes(ymin = cal_lwr, ymax = cal_upr, fill = "Calibrated"), alpha = 0.15) +
  geom_ribbon(aes(ymin = err_lwr, ymax = err_upr, fill = "Error-injected"), alpha = 0.20) +
  geom_line(aes(y = brier_cal, colour = "Calibrated")) +
  geom_line(aes(y = brier_err, colour = "Error-injected")) +
  geom_point(aes(y = brier_cal, colour = "Calibrated")) +
  geom_point(aes(y = brier_err, colour = "Error-injected")) +
  scale_colour_manual(name = NULL,
                      values = c("Calibrated" = "darkred", "Error-injected" = "steelblue")) +
  scale_fill_manual(name = NULL,
                    values = c("Calibrated" = "darkred", "Error-injected" = "steelblue")) +
  guides(fill = "none") +
  labs(x = "Base rate", y = "Brier score") +
  theme_bw() +
  theme(legend.position = "top")
print(brier_effect_plot2)

# Calculate and plot the delta
boot_brier$gap <- boot_brier$brier_err - boot_brier$brier_cal
brier_delta_plot <- ggplot(boot_brier, aes(base_rate, gap)) +
  geom_line() +
  geom_point() +
  labs(x = "Base rate", y = "Brier inflation from error",
       title = "Delta: Calibrated vs. Error-injected") +
  theme_bw()
print(brier_delta_plot)
ggarrange(brier_effect_plot2, brier_delta_plot)

## END ##