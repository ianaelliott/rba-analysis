#------------------------------------------------------------------------------#
## HOLLIS RISK ANALYSIS EXPLORATION v11 ##
#------------------------------------------------------------------------------#

library(epitools)
library(furrr)
library(pwr)
library(effectsize)
library(tidyr)
library(data.table)
library(future.apply)
library(future)
library(ggpubr)
library(tidyverse)

#rm(list=ls()) # Clear environment
options(scipen = 999) # Disable scientific notation

seed1 <- 6012021  # Alfie
seed2 <- 3222023  # Beya
seed3 <- 9061991  # Bizzle

#------------------------------------------------------------------------------#
## 1. VALIDATE THE HOLLIS ANALYSIS BEFORE REPLICATION ##
#------------------------------------------------------------------------------#

# Parameters
model_auc <- 0.813       # ogrs v2 auc used by Hollis
hollis_n  <- 12924       # original study sample size
hollis_obs_rx <- 0.619   # observed recidivism rate
hollis_exp_rx <- 0.696   # expected counterfactual rate

hollis_obs <- round(hollis_n*hollis_obs_rx,0)
hollis_pred <- hollis_n*hollis_exp_rx

# Contingency table
hollis_contingency <- matrix (
  c(hollis_n - hollis_obs,   # observed nrx
    hollis_n - hollis_pred,  # predicted nrx
    hollis_obs,              # observed rx
    hollis_pred              # predicted rx
  ),
  nrow = 2,
  dimnames = list(
    c("Observed", "Predicted"),
    c("Non-recidivist", "Recidivist")
  )
)
print(hollis_contingency)

# Statistical tests 
chi_hollis <- chisq.test(hollis_contingency)

# Goodness-of-fit and risk ratio
chi_hollis_gof <- effectsize::cohens_w(chi_hollis, adjust = FALSE)[,1]
chi_hollis_rr <- epitools::riskratio(hollis_contingency)

# View outcomes
cat("Contingency table outcomes:\n",
    "- Chi-squared:    ", round(chi_hollis$statistic, 3), "\n",
    "- Chi-squared df: ", round(chi_hollis$parameter, 3), "\n",
    "- Chi-squared p:  ", round(chi_hollis$p.value, 3), "\n",
    "- Risk ratio:     ", round(chi_hollis_rr$measure[2,1], 3), "\n",
    "- RR 95% CI:      [", round(chi_hollis_rr$measure[2,2], 3),", ",round(chi_hollis_rr$measure[2,3], 3),"]", "\n",
    "- Cohen's w:      ", round(chi_hollis_gof, 3),
    sep = "")

#------------------------------------------------------------------------------#
## 2. REPRODUCIBLE SYNTHETIC DATASET WITH PLAUSIBLE PARAMETERS ##
#------------------------------------------------------------------------------#
# Generate beta-distributed risk scores

# Parameters tuned for approximate realism (mean ~0.7, sd ~0.14-0.18)
scaling_factor <- 10  # lower = more dispersion

set.seed(seed1)
synth_hollis <- data.frame(
  risk = rbeta(
    n = hollis_n,
    shape1 = 0.696 * scaling_factor,
    shape2 = (1 - 0.696) * scaling_factor
  )
) 

# Validate parameters approximate Hollis's
cat(
  "Approximate Mean:", round(mean(synth_hollis$risk), 3), # 0.695
  "\nApproximate SD  :", round(sd(synth_hollis$risk), 3)  # 0.139
)

# Visualize distribution
hist(synth_hollis$risk, 
     main = "Synthetic Risk Scores", 
     xlab = "Predicted Recidivism Probability")

# Write data for future reference
saveRDS(synth_hollis, "synth_hollis.rds")

#------------------------------------------------------------------------------#
## 3. NAIVE RISK SCORE ANALYSIS (RSA) SIMULATION ##
#------------------------------------------------------------------------------#

# Generate simulated "true" outcomes
set.seed(seed2)
synth_outcomes <- synth_hollis %>% 
  dplyr::mutate(
    rx = rbinom(n(), 1, risk)  # "true" outcomes from risk scores
  )

# Essential metrics calculation
hollis_counterfactual <- hollis_contingency[2,2]    # 8995
simulated_counterfactual <- sum(synth_outcomes$rx)  # 8991

# Validation output
metrics <- list(
  # core counterfactual comparison
  absolute_bias = hollis_counterfactual - simulated_counterfactual,
  relative_bias = (hollis_counterfactual - simulated_counterfactual) / hollis_counterfactual,
  # distribution validation
  risk_mean = mean(synth_outcomes$risk),
  risk_sd = sd(synth_outcomes$risk),
  # reconviction rates
  simulated_rate = simulated_counterfactual / hollis_n,
  hollis_rate = 0.696
)

cat(
  "Synthetic data validation\n",
  "Risk Score Fidelity:\n",
  "- Target Mean: 0.696 | Achieved: ", round(metrics$risk_mean, 3), "\n",
  "- Target SD: ~0.13-0.18 | Achieved: ", round(metrics$risk_sd, 3), "\n\n",
  "Reconviction Comparison:\n",
  "- Simulated: ", simulated_counterfactual, " (", round(metrics$simulated_rate*100,1), "%)\n",
  "- Hollis:    ", hollis_counterfactual, " (69.6%)\n",
  "- Delta:     ", metrics$absolute_bias, 
  " (", round(metrics$relative_bias*100,3), "% relative bias)\n",
  sep = ""
)

# Statistical validation (comparison to Hollis's reported rate)
rate_test <- binom.test(
  x = simulated_counterfactual,
  n = hollis_n,
  p = 0.696
)

cat(
  "\nBinomial Test vs. Hollis Rate:\n",
  "- p-value: ", signif(rate_test$p.value, 3), "\n",
  "- 95% CI: [", 
  round(rate_test$conf.int[1], 3), ", ", 
  round(rate_test$conf.int[2], 3), "]\n",
  sep = ""
)
# Hollis' rate (0.696) falls within the simulation's CI 
# Synthetic data plausibly matches original parameters

#--------------------------------#
## SENSITIVITY CHECKS

# 1. Ensure changes in SD != large change in predicted outcomes
risk_sd_ci <- quantile(synth_hollis$risk, probs = c(0.025, 0.975)) %>% 
  round(3)

cat(
  "Risk Score SD Validation:\n",
  "- Simulated SD: ", round(sd(synth_hollis$risk), 3), 
  " 95% CI [", risk_sd_ci[1], ", ", risk_sd_ci[2], "]\n", # wide CI reflects beta distribution's skew
  "- Brinn et al. SD: 0.18\n",
  sep = ""
)

# Test higher SD (Brinn's 0.18) by adjusting scaling factor 
set.seed(seed3)
synth_hollis_high_sd <- data.frame(
  risk = rbeta(12924, 0.696*5, (1 - 0.696)*5)  #scaling factor = 5; SD ~0.18
)

synth_outcomes_high_sd <- synth_hollis_high_sd %>% 
  dplyr::mutate(
    rx = rbinom(n(), 1, risk)  # true outcomes from risk scores
  )

# Compare reconvictions between SD = 0.13 vs. SD = 0.18
cat(
  "Hollis's counterfactual:  ", hollis_counterfactual,            # 8970
  "\nSimulated counterfactual: ", simulated_counterfactual,       # 8991
  "\nHigh SD counterfactual:   ", sum(synth_outcomes_high_sd$rx), # 8971
  "\nRelative over-prediction: " , round((simulated_counterfactual-hollis_counterfactual)/hollis_counterfactual,5)
)

# 2. Calculate relative over-prediction for comparison-to-observed fairness
overprediction_table <- tibble(
  Study = c("Hollis", "Synthetic Hollis"),
  Predicted = c(8995, 8991),
  Observed = rep(8000, 2),
  Overprediction = Predicted - Observed,
  Overprediction_pc = (Predicted / Observed - 1) * 100
)
print(overprediction_table)

# Proportional difference between studies
overprediction_ratio <- overprediction_table$Overprediction_pc[2] / overprediction_table$Overprediction_pc[1]  # ~0.996

# Two proportion z-test for divergence (Hollis's vs our prediction)
prop_test <- prop.test(
  x = c(simulated_counterfactual, hollis_counterfactual),
  n = c(hollis_n, hollis_n)
)
cat("z =", round(prop_test$statistic, 2), ", p =", prop_test$p.value)

# Test +/-1% variation in observed counts
overprediction_sensitivity <- map_dfr(
  c(7920, 8000, 8080),  # -1%, base, +1%
  \(obs) {
    tibble::tibble(
      Observed = obs,
      Overprediction = 8991 - obs,
      Overprediction_pc = (8991 / obs - 1) * 100
    )
  }
)
print(overprediction_sensitivity)

# Save for reproducibility
saveRDS(
  list(
    data = synth_outcomes,
    metrics = metrics,
    test = rate_test
  ),
  "synth_outcomes.rds"
)

#------------------------------------------------------------------------------#
## 4. MEASUREMENT ERROR CORRECTION ##
#------------------------------------------------------------------------------#

# Apply AUC-based measurement error correction
set.seed(seed3) 
synth_error <- synth_outcomes %>% 
  dplyr::mutate(
    rx_auc = ifelse(runif(n()) < model_auc, rx, 1 - rx) # vectorized error injection
  )

# Validate error injection
error_matrix <- synth_error %>%
  count(rx, rx_auc) %>% 
  dplyr::mutate(
    category = case_when(
      rx == 1 & rx_auc == 1 ~ "True Positive",
      rx == 0 & rx_auc == 0 ~ "True Negative",
      rx == 1 & rx_auc == 0 ~ "False Negative",
      rx == 0 & rx_auc == 1 ~ "False Positive"
    )
  ) %>% 
  dplyr::group_by(category) %>% 
  dplyr::summarise(n = sum(n), .groups = "drop") %>% 
  dplyr::mutate(prop = n/nrow(synth_error))
print(error_matrix)

# Isolate total proportion flipped
error_matrix %>%
  dplyr::filter(category == "False Positive" | category == "False Negative") %>%
  dplyr::summarise(sum = sum(prop))  # 0.189 (18.9.1%)

# Isolate total proportion maintained
error_matrix %>%
  dplyr::filter(category == "True Positive" | category == "True Negative") %>%
  dplyr::summarise(sum = sum(prop))  # 0.811 (81.1%)

# Contingency table
hollis_error_contingency <- matrix (
  c(hollis_n - hollis_obs,               # observed nrx
    hollis_n - sum(synth_error$rx_auc),  # adjusted nrx
    hollis_obs,                          # observed rx
    sum(synth_error$rx_auc)              # adjusted rx
  ),
  nrow = 2,
  dimnames = list(
    c("Observed", "Predicted"),
    c("Non-recidivist", "Recidivist")
  )
)
print(hollis_error_contingency)  # 8093

# Chi-square & Cohen’s w for effect size
chi_adjusted <- chisq.test(hollis_error_contingency)  

# Goodness-of-fit and risk ratio
chi_adjusted_gof <- round(effectsize::cohens_w(chi_adjusted, adjust = FALSE)[,1],3)
chi_adjusted_rr <- epitools::riskratio(hollis_error_contingency) 

# View outcomes
cat("Contingency table outcomes:\n",
    "- Chi-squared:    ", round(chi_adjusted$statistic, 3), "\n",
    "- Chi-squared df: ", round(chi_adjusted$parameter, 3), "\n",
    "- Chi-squared p:  ", round(chi_adjusted$p.value, 3), "\n",
    "- Risk ratio:     ", round(chi_adjusted_rr$measure[2,1], 3), "\n",
    "- RR 95% CI:      [", round(chi_adjusted_rr$measure[2,2], 3),", ",round(chi_adjusted_rr$measure[2,3], 3),"]", "\n",
    "- Cohen's w:      ", round(chi_adjusted_gof, 3),
    sep = "")

# Save for reproducibility
saveRDS(synth_error, "synth_error.rds")

#------------------------------------------------------------------------------#
## 5. BOOTSTRAP SIMULATION OF NAIVE vs. CORRECTED OUTCOMES ##
#------------------------------------------------------------------------------#

# Simulation parameters
n_boot <- 50000     # number of bootstrap iterations
confidence <- 0.95  # confidence level for intervals

# Simulation function
simulate_prediction <- function(n_boot, df, auc = NULL) {
    future_replicate(
    n_boot,
    expr = {
      if (is.null(auc)) {
        # naive RSA: generate outcomes from risk scores
        sum(rbinom(nrow(df), 1, df$risk))
      } else {
        # measurement error correction
        flips <- runif(nrow(df)) > auc
        sum(ifelse(flips, 1 - df$rx, df$rx))
      }
    },
    future.seed = TRUE,   # critical for reproducibility
    future.stdout = FALSE # cleaner output
  )
}

# Run simulations
plan(multisession, workers = availableCores() - 2) # preserve 2 cores for system

set.seed(seed1) 
system.time({
  naive_boot <- simulate_prediction(n_boot, synth_outcomes) # runs for naive outcomes
  corrected_boot <- simulate_prediction(n_boot, synth_outcomes, model_auc) # runs for adjusted outcomes
})

# Prepare data for visualization
boot_outcomes <- data.frame(
  Method = rep(c("Naive", "Adjusted"), each = length(naive_boot)),
  Count = c(naive_boot, corrected_boot)
)

# Save for reproducibility
saveRDS(boot_outcomes, "boot_outcomes.rds")

# Calculate summary stats for annotations
boot_summary <- boot_outcomes %>%
  dplyr::group_by(Method) %>%
  dplyr::summarise(
    mean = mean(Count),
    ci_low = quantile(Count, 0.025),
    ci_high = quantile(Count, 0.975)
  )

# Visualize the bootstrapped data
boot_plot <- ggplot(boot_outcomes, aes(x = Count, fill = Method)) +
  geom_density(alpha = 0.6, adjust = 1.5) +
  geom_rect(data = boot_summary, 
            aes(xmin = ci_low, xmax = ci_high, ymin = 0, ymax = Inf, fill = Method), 
            alpha = 0.15, 
            inherit.aes = FALSE) +
  geom_vline(data = boot_summary, 
             aes(xintercept = mean, color = Method),
             linetype = "dashed", linewidth = 0.6,
             alpha = 0.9) +
  scale_fill_viridis_d(begin = 0.2, end = 0.5) +
  scale_color_viridis_d(begin = 0.2, end = 0.5) +
  labs(
    x = "Number of reoffenders (RSA)", 
    y = "Density") +
  #xlim(7500, 9500) +
  theme_bw() +
  theme(legend.position = "top")
print(boot_plot)

# Consolidated statistical analyses
naive_mean <- mean(naive_boot)
adjusted_mean <- mean(corrected_boot)

boot_results <- list(
  naive = list(
    mean = naive_mean,
    sd = sd(naive_boot),
    ci = quantile(naive_boot, c(0.025, 0.975))
  ),
  corrected = list(
    mean = adjusted_mean,
    sd = sd(corrected_boot),
    ci = quantile(corrected_boot, c(0.025, 0.975))
  ),
  comparison = list(
    delta_abs = mean(naive_boot - corrected_boot),
    delta_rel = mean((naive_boot - corrected_boot)/naive_boot),
    
    contingency = matrix(c(
      naive_mean, hollis_n - naive_mean,
      adjusted_mean, hollis_n - adjusted_mean
    ), nrow = 2),
    
    risk_ratio = epitools::riskratio(matrix(c(
      naive_mean, hollis_n - naive_mean,
      adjusted_mean, hollis_n - adjusted_mean
    ), nrow = 2)),
    
    cohens_w = effectsize::cohens_w(
      matrix(c(
        naive_mean, hollis_n - naive_mean,
        adjusted_mean, hollis_n - adjusted_mean
      ), nrow = 2)
    )
  )
)

# Results output
cat("Key Findings RSA:\n",
    "Naive Model (error-free assumption):\n",
    "  Mean = ", round(boot_results$naive$mean), "\n",
    "  SD = ", round(boot_results$naive$sd, 1), "\n",
    "  95% CI = [", boot_results$naive$ci[1], ", ", boot_results$naive$ci[2], "]\n\n",
    
    "Corrected Model (AUC = ", model_auc, "):\n",
    "  Mean = ", round(boot_results$corrected$mean), "\n",
    "  SD = ", round(boot_results$corrected$sd, 1), "\n",
    "  95% CI = [", boot_results$corrected$ci[1], ", ", boot_results$corrected$ci[2], "]\n\n",
    
    "Comparison Metrics:\n",
    "  Absolute delta: ", round(boot_results$comparison$delta_abs), "\n",
    "   - [95% CI: ", round(quantile(naive_boot - corrected_boot, 0.025)), 
    ", ", round(quantile(naive_boot - corrected_boot, 0.975)), "]\n",
    "  Relative delta: ", round(boot_results$comparison$delta_rel*100, 1), "%\n",
    "  Risk Ratio: ", round(boot_results$comparison$risk_ratio$measure[2,1], 2), 
    " (95% CI: ", round(boot_results$comparison$risk_ratio$measure[2,2], 2), "-",
    round(boot_results$comparison$risk_ratio$measure[2,3], 2), ")\n",
    "  Cohen's w: ", round(boot_results$comparison$cohens_w[[1]], 2), 
    sep = "")

## How often would we mistakenly detect a "significant" effect if we ignore 
## measurement error? Does correction reduce false-positive rates?

# For naive bootstrap
naive_pvals <- map_dbl(naive_boot, ~ {
  cont_table <- matrix(c(.x, hollis_n - .x, 8000, hollis_n - 8000), nrow = 2)
  chisq.test(cont_table)$p.value
})
mean(naive_pvals < 0.05)  # 1.0 (100%)
min(naive_pvals); max(naive_pvals)
table(naive_pvals < 0.05)

# For corrected bootstrap
corrected_pvals <- map_dbl(corrected_boot, ~ {
  cont_table <- matrix(c(.x, hollis_n - .x, 8000, hollis_n - 8000), nrow = 2)
  chisq.test(cont_table)$p.value
})
mean(corrected_pvals < 0.05)  # 0.0716 (7.2%)
min(corrected_pvals); max(corrected_pvals)
table(corrected_pvals < 0.05)

# Visualized
all_pvals <- data.frame(
  Naive = naive_pvals,
  Adjusted = corrected_pvals 
)

all_pvals %>%
  tidyr::pivot_longer(cols = c(Naive, Adjusted), 
               names_to = "Analysis",
               values_to = "p") %>%
  ggplot(., aes(x = p, fill = Analysis)) +
  geom_histogram(position="dodge") +
  scale_fill_viridis_d(begin = 0.2, end = 0.5) +
  theme_bw() +
  theme(legend.position = "top")

#------------------------------------------------------------------------------#
## 6. SENSITIVITY ANALYSIS: SYSTEMATIC PARAMETER VARIATION ##
#------------------------------------------------------------------------------#

# Sensitivity parameters 
baseline_levels <- c(0.100, 0.250, 0.500, 0.695)  # low, medium, median, and higher 
auc_levels <- c(0.713, 0.813, 0.913)  # poor, typical, excellent

# Sensitivity analyses simulation function
simulate_sensitivities <- function(auc, baseline, hollis_n = 12924, n_boot = 50000) {
  
  # flexible beta parameterization
  get_beta_shape <- function(mu, scaling = 10) {
    list(shape1 = mu * scaling, shape2 = (1 - mu) * scaling)
  }
  params <- get_beta_shape(baseline)
  risk_scores <- rbeta(hollis_n, params$shape1, params$shape2)
  
  # vectorized outcome generation
  true_outcomes <- rbinom(hollis_n, 1, risk_scores)
  flips <- runif(hollis_n) > auc
  corrected_outcomes <- ifelse(flips, 1 - true_outcomes, true_outcomes)
  
  # efficient bootstrap
  future_replicate(
    n_boot,
    expr = {
      idx <- sample(hollis_n, replace = TRUE)
      c(naive = sum(true_outcomes[idx]), corrected = sum(corrected_outcomes[idx]))
    },
    future.seed = TRUE
  ) %>% t() %>% as_tibble()
}

# Run sensitivity simulations
plan(multisession, workers = availableCores() - 2)

set.seed(seed2) 
system.time({
  sensitivity_outcomes <- expand.grid(auc = auc_levels, baseline = baseline_levels) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      simulations = future_pmap(
        list(auc, baseline),
        simulate_sensitivities,
        .options = furrr_options(seed = TRUE, stdout = FALSE),
        .progress = TRUE
      )
    )
})

# Save for reproducibility
saveRDS(sensitivity_outcomes, "sensitivity_outcomes.rds")

# Consolidated statistical analyses
sensitivity_results <- sensitivity_outcomes %>% 
  tidyr::unnest(simulations) %>% 
  dplyr::group_by(auc, baseline) %>% 
  dplyr::summarise(
    across(c(naive, corrected), 
           list(mean = mean, sd = sd, 
                ci_low = ~quantile(.x, 0.025), 
                ci_high = ~quantile(.x, 0.975))),
    .groups = "drop"
  ) %>% 
  dplyr::mutate(
    absolute_bias = naive_mean - corrected_mean,
    relative_bias = absolute_bias / naive_mean,
    risk_ratio = naive_mean / corrected_mean,
    cohens_w = map2_dbl(
      naive_mean, corrected_mean,
      ~effectsize::cohens_w(
        matrix(c(
          .x, 12924 - .x,  # using the fixed RBA sample size
          .y, 12924 - .y
        ), nrow = 2)
      )[[1]]
    )
  )

# Print output
sensitivity_results %>% 
  tibble::as_tibble() %>% 
  print(n = Inf, width = Inf)

# Save results
write.csv(sensitivity_results, "sensitivity_results.csv") 

# Effect size by AUC and baseline
effect_plot_rsa <- ggplot(sensitivity_results, 
                          aes(factor(auc), risk_ratio, color = factor(baseline), group = factor(baseline))) +
  geom_point(size = 2) + 
  geom_line() +
  geom_errorbar(aes(
    ymin = naive_ci_low/corrected_ci_high, 
    ymax = naive_ci_high/corrected_ci_low),
    width = 0.05) + 
  ylim(0, 1.25) +
  scale_color_viridis_d(
    name = "Baseline rate", 
    option = "viridis", 
    begin = 0, 
    end = 0.9) +
  labs(
    x = "AUC", 
    y = "Risk ratio (naive/corrected)") +
  theme_bw() + 
  theme(legend.position = "top")
print(effect_plot_rsa)

# Bias distribution plot
bias_plot_rsa <- ggplot(sensitivity_results, 
                        aes(factor(baseline), relative_bias, fill = factor(auc))) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d(name = "AUC",
                       option = "viridis", 
                       begin = 0, 
                       end = 0.9) +
  labs(
    x = "Baseline rate", 
    y = "Relative bias") +
  ylim(-3.5, 0.5) +
  theme_bw() + 
  theme(legend.position = "top")
print(bias_plot_rsa)

# Combine with RSA plots
ggarrange(effect_plot_rsa, bias_plot_rsa)

#------------------------------------------------------------------------------#
## 7. RBA SENSITIVITY ANALYSES: INDIVIDUAL-LEVEL COARSENING ##
#------------------------------------------------------------------------------#
# We'll go straight to the sensitivity analyses here as they will include our
# original 0.695 and 0.813 parameters. We proportionally adjust the original 
# reference rates (for 0.695) to other baselines, preserving between-band ratios

# Calculate the band rates for AUC = 0.813 and baseline 0.695 via synth_hollis.
# This uses published cut-offs for OGRS3 (not OGRS2) based on a similar population. 
original_rates <- synth_hollis %>%
  dplyr::mutate(
    band = cut(
      risk,
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("Low", "Medium", "High", "Very High"), 
      include.lowest = TRUE)) %>%
  dplyr::group_by(band) %>%
  dplyr::summarise(mean = mean(risk)) %>%
  deframe()

# Function to scale rates proportionally
scale_rates <- function(new_baseline, origin_rates, origin_mean) {
  ratio <- origin_mean / origin_rates
  scaled_rates <- new_baseline / ratio
  pmin(pmax(scaled_rates, 0), 1)
}

# Generate the exogenous rates
exogenous_rates <- sapply(baseline_levels, function(x) {
  scale_rates(x, original_rates, 0.695)
})

# Generate a matrix of plausible rates
exogenous_rates_matrix <- matrix(
  exogenous_rates,
  nrow = 4,
  dimnames = list(
    names(original_rates),
    paste0("baseline_", sprintf("%.3f", baseline_levels))
  )
)
print(exogenous_rates_matrix)

# Function to scale thresholds proportionally
scale_breaks <- function(new_baseline, origin_breaks, origin_mean) {
  ratio <- origin_mean / origin_breaks
  scaled_rates <- new_baseline / ratio
  pmin(pmax(scaled_rates, 0), 1)
}

# Specify the original OGRS3 breaks (assuming also approximately true for OGRS2)
original_breaks <- c(0, 0.25, 0.5, 0.75, 1)

# Generate a matrix of plausible breaks
exogenous_breaks <- sapply(baseline_levels, function(x) {
  scaled_breaks <- original_breaks * (x / 0.695)
  clamped <- pmin(pmax(scaled_breaks, 0), 1)
  clamped[5] <- 1  # ensure final break is always 1
  sort(clamped)    # force ascending order
})

colnames(exogenous_breaks) <- paste0("baseline_", sprintf("%.3f", baseline_levels))
rownames(exogenous_breaks) <- c("Break0", "Break1", "Break2", "Break3", "Break4")
print(exogenous_breaks)

# Sensitivity analyses RBA simulation function
simulate_rba_sensitivities <- function(auc, baseline, hollis_n = 12924, n_boot = 50000) {
  # generate consistent baseline label
  baseline_label <- paste0("baseline_", sprintf("%.3f", baseline))
  
  # parameter retrieval
  breaks <- exogenous_breaks[, baseline_label]
  rates <- exogenous_rates_matrix[, baseline_label]
  rates <- pmax(pmin(rates, 1), 0)  # force within [0,1]
  
  # latent risk generation
  get_beta_shape <- function(mu, scaling = 10) {
    list(shape1 = mu * scaling, shape2 = (1 - mu) * scaling)
  }
  
  params <- get_beta_shape(baseline)
  risk_scores <- rbeta(hollis_n, params$shape1, params$shape2)
  
  # risk band assignment
  band_assignments <- cut(
    risk_scores,
    breaks = breaks,
    labels = names(rates),
    include.lowest = TRUE
  ) %>% 
    as.character() %>%  # convert factor to character for safer handling
    
    # replace any NA assignments with lowest risk category
    replace_na(names(rates)[1]) %>% 
    factor(levels = names(rates))
  
  # rate assignment validation
  coarsened_rates <- rates[band_assignments]
  
  # outcome generation
  true_outcomes <- rbinom(hollis_n, 1, coarsened_rates)
  
  # measurement error injection
  flips <- runif(hollis_n) > auc
  corrected_outcomes <- ifelse(flips, 1 - true_outcomes, true_outcomes)
  
  # bootstrap
  future_replicate(
    n_boot,
    expr = {
      idx <- sample(hollis_n, replace = TRUE)
      c(naive = sum(true_outcomes[idx]), corrected = sum(corrected_outcomes[idx]))
    },
    future.seed = TRUE
  ) %>% t() %>% as_tibble()
}

# Parallelized execution
plan(multisession, workers = availableCores() - 2)
set.seed(seed2)

system.time({
  sensitivity_outcomes_rba <- expand.grid(
    auc = auc_levels,
    baseline = baseline_levels # Use pre-formatted levels
  ) %>%
    dplyr::mutate(
      simulations = future_pmap(
        list(auc, baseline),
        simulate_rba_sensitivities,
        .options = furrr_options(seed = TRUE, stdout = FALSE),
        .progress = TRUE
      )
    )
})

# Save for reproducibility
saveRDS(sensitivity_outcomes_rba, "sensitivity_outcomes_rba.rds")

sensitivity_results_rba <- sensitivity_outcomes_rba %>% 
  tidyr::unnest(simulations) %>% 
  dplyr::group_by(auc, baseline) %>% 
  dplyr::summarise(
    across(c(naive, corrected), 
           list(mean = mean, sd = sd, 
                ci_low = ~quantile(.x, 0.025), 
                ci_high = ~quantile(.x, 0.975))),
    .groups = "drop"
  ) %>% 
  dplyr::mutate(
    absolute_bias = naive_mean - corrected_mean,
    relative_bias = absolute_bias / naive_mean,
    risk_ratio = naive_mean / corrected_mean,
    cohens_w = map2_dbl(
      naive_mean, corrected_mean,
      ~effectsize::cohens_w(
        matrix(c(
          .x, 12924 - .x,
          .y, 12924 - .y
        ), nrow = 2)
      )[[1]]
    )
  )

# Print formatted results
sensitivity_results_rba %>% 
  dplyr::mutate(across(where(is.numeric), ~round(., 3))) %>% 
  print(n = Inf, width = Inf)

# Save results
write.csv(sensitivity_results_rba, "sensitivity_results_rba.csv")

# Effect size plot for RBA IC
effect_plot_rba <- ggplot(sensitivity_results_rba, 
                          aes(factor(auc), risk_ratio, 
                              color = factor(baseline), 
                              group = factor(baseline))) +
  geom_point(size = 2) + 
  geom_line() +
  geom_errorbar(aes(
    ymin = naive_ci_low/corrected_ci_high,
    ymax = naive_ci_high/corrected_ci_low),
    width = 0.05) + 
  ylim(0, 1.25) +
  scale_color_viridis_d(
    name = "Baseline recidivism rate", 
    option = "viridis",
    begin = 0, 
    end = 0.9) +
  labs(
    x = "AUC", 
    y = "Risk ratio (naive/corrected)",
    title = "Risk band analysis (OL)") +
  theme_bw() +
  theme(legend.position = "top")
print(effect_plot_rba)

# Bias plot for RBA IC
bias_plot_rba <- ggplot(sensitivity_results_rba, 
                        aes(factor(baseline), absolute_bias, 
                            fill = factor(auc))) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d(
    name = "AUC",
    option = "viridis",
    begin = 0,
    end = 0.9) +
  labs(
    x = "Baseline recidivism rate", 
    y = "Absolute bias (naive - corrected)",
    #title = "Risk band analysis (OL)"
    ) +
  theme_bw() +
  theme(legend.position = "top")
print(bias_plot_rba)

# Combine with RSA plots
ggarrange(effect_plot_rsa, effect_plot_rba, ncol = 2, common.legend = TRUE)
ggarrange(bias_plot_rsa, bias_plot_rba, ncol = 2, common.legend = TRUE)

# Plot the difference between naive and corrected mean counterfactuals
ggplot(sensitivity_results_rba , aes(x = factor(baseline), y = absolute_bias , color = factor(auc), group = factor(auc))) +
  geom_point(size = 2) +
  geom_line() +
  scale_color_viridis_d(name = "AUC", 
                        option = "viridis", 
                        begin = 0, 
                        end = 0.9) +
  labs(x = "Baseline recidivism rate", 
       y = "Difference (naive - corrected)", 
       title = "Impact of coarsening on mean counterfactual") +
  theme_bw() +
  theme(legend.position = "top")

# Bias distribution plot
bias_plot_rba_ic <- ggplot(sensitivity_results_rba, 
                        aes(factor(baseline), relative_bias, fill = factor(auc))) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d(name = "AUC",
                       option = "viridis", 
                       begin = 0, 
                       end = 0.9) +
  labs(
    x = "Baseline rate", 
    y = "Relative bias") +
  ylim(-3.5, 0.5) +
  theme_bw() + 
  theme(legend.position = "top")
print(bias_plot_rba_ic)

ggarrange(bias_plot_rsa, bias_plot_rba_ic, ncol = 2, common.legend = TRUE)

#------------------------------------------------------------------------------#
## 8. RBA SENSITIVITY ANALYSES: CLASSIFICATION-LEVEL COARSENING ##
#------------------------------------------------------------------------------#

simulate_rba_classification <- function(auc, baseline, hollis_n = 12924, n_boot = 50000) {
  # Get parameters for current baseline
  baseline_label <- paste0("baseline_", sprintf("%.3f", baseline))
  breaks <- exogenous_breaks[, baseline_label]
  rates <- setNames(exogenous_rates_matrix[, baseline_label], rownames(exogenous_rates_matrix))
  
  # Generate latent risk scores
  get_beta_shape <- function(mu, scaling = 10) {
    list(shape1 = mu * scaling, shape2 = (1 - mu) * scaling)
  }
  params <- get_beta_shape(baseline)
  risk_scores <- rbeta(hollis_n, params$shape1, params$shape2)
  
  # Assign original bands
  original_bands <- cut(
    risk_scores,
    breaks = breaks,
    labels = names(rates),
    include.lowest = TRUE
  )
  
  # Generate flipped classifications
  flipped_bands <- original_bands
  flip_mask <- runif(hollis_n) > auc
  flipped_bands[flip_mask] <- sample(
    names(rates), 
    sum(flip_mask), 
    replace = TRUE
  )
  
  # Create frequency tables
  naive_table <- table(original_bands)
  corrected_table <- table(flipped_bands)
  
  # Ensure all bands are represented (in case of empty categories)
  full_bands <- factor(levels(original_bands), levels = levels(original_bands))
  naive_table <- table(factor(original_bands, levels = levels(full_bands)))
  corrected_table <- table(factor(flipped_bands, levels = levels(full_bands)))
  
  # Calculate predicted reoffenders
  naive_predicted <- sum(naive_table * rates)
  corrected_predicted <- sum(corrected_table * rates)
  
  # Bootstrap resampling
  future_replicate(
    n_boot,
    expr = {
      # Resample with replacement
      idx <- sample(hollis_n, replace = TRUE)
      
      # Create bootstrapped tables
      bs_naive <- table(factor(original_bands[idx], levels = levels(full_bands)))
      bs_corrected <- table(factor(flipped_bands[idx], levels = levels(full_bands)))
      
      c(
        naive = sum(bs_naive * rates),
        corrected = sum(bs_corrected * rates)
      )
    },
    future.seed = TRUE
  ) %>% t() %>% as_tibble()
}

# Execute parallel simulation
plan(multisession, workers = availableCores() - 2)
set.seed(seed2)

system.time({
  sensitivity_rba_class <- expand.grid(
    auc = auc_levels,
    baseline = baseline_levels
  ) %>%
    dplyr::mutate(
      simulations = future_pmap(
        list(auc, baseline),
        simulate_rba_classification,
        .options = furrr_options(seed = TRUE, stdout = FALSE),
        .progress = TRUE
      )
    )
})

# Save for reproducibility
saveRDS(sensitivity_rba_class, "sensitivity_rba_class.rds")

# Analysis and visualization
sensitivity_results_rba_class <- sensitivity_rba_class %>% 
  unnest(simulations) %>% 
  dplyr::group_by(auc, baseline) %>% 
  dplyr::summarise(
    across(c(naive, corrected), 
           list(mean = ~mean(.x, na.rm = TRUE), 
                sd = ~sd(.x, na.rm = TRUE), 
                ci_low = ~quantile(.x, 0.025, na.rm = TRUE), 
                ci_high = ~quantile(.x, 0.975, na.rm = TRUE))),
    .groups = "drop"
  ) %>% 
  dplyr::mutate(
    absolute_bias = naive_mean - corrected_mean,
    relative_bias = absolute_bias / naive_mean,
    risk_ratio = naive_mean / corrected_mean,
    cohens_w = map2_dbl(
      naive_mean, corrected_mean,
      ~effectsize::cohens_w(
        matrix(c(
          .x, 12924 - .x,
          .y, 12924 - .y
        ), nrow = 2)
      )[[1]]
    )
  )

# Print formatted results
sensitivity_results_rba_class %>% 
  dplyr::mutate(across(where(is.numeric), ~round(., 3))) %>% 
  print(n = Inf, width = Inf)

# Save results to a file for transfer to MS Word
write.csv(sensitivity_results_rba_class, "sensitivity_results_rba_class.csv")

# Generate comparative plots
effect_plot_rba_class <- ggplot(sensitivity_results_rba_class, 
                                aes(factor(auc), risk_ratio, 
                                    color = factor(baseline), 
                                    group = factor(baseline))) +
  geom_point(size = 2) + 
  geom_line() +
  geom_errorbar(aes(
    ymin = naive_ci_low/corrected_ci_high,
    ymax = naive_ci_high/corrected_ci_low),
    width = 0.05) + 
  ylim(0, 1.25) +
  scale_color_viridis_d(
    name = "Baseline recidivism rate", 
    option = "viridis",
    begin = 0, 
    end = 0.9) +
  labs(
    x = "AUC", 
    y = "Risk ratio (naive/corrected)",
    title = "Risk band analysis (CL)") +
  theme_bw() +
  theme(legend.position = "top")
print(effect_plot_rba_class)

# Bias distribution plot
bias_plot_rba_ac <- ggplot(sensitivity_results_rba_class, 
                           aes(factor(baseline), relative_bias, fill = factor(auc))) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d(name = "AUC",
                       option = "viridis", 
                       begin = 0, 
                       end = 0.9) +
  labs(
    x = "Baseline rate", 
    y = "Relative bias") +
  ylim(-3.5, 0.5) +
  theme_bw() + 
  theme(legend.position = "top")
print(bias_plot_rba_ac)

# Combine all plots
ggarrange(effect_plot_rsa, effect_plot_rba, effect_plot_rba_class, 
          ncol = 3, common.legend = TRUE)
ggarrange(bias_plot_rsa, bias_plot_rba_ic, bias_plot_rba_ac, 
          ncol = 3, common.legend = TRUE)

#------------------------------------------------------------------------------#
## 9. COMBINED SENSITIVITY ANALYSES EXPLORATION ##
#------------------------------------------------------------------------------#

sensitivity_results
sensitivity_results_rba
sensitivity_results_rba_class

all_sensitivities <- bind_rows(
  list(
    RSA = sensitivity_results,
    RBA_OL = sensitivity_results_rba,
    RBA_CL = sensitivity_results_rba_class
  ),
  .id = "analysis_type"
)

all_sensitivities %>%
  dplyr::select(analysis_type, naive_mean, corrected_mean, baseline, auc) %>%
  tidyr::pivot_longer(
    cols = c(naive_mean, corrected_mean),  # Explicitly define columns to pivot
    names_to = "injection",
    values_to = "mean"
  ) %>%
  dplyr::mutate(analysis_type = factor(analysis_type, levels = c("RSA","RBA_OL","RBA_CL")),
         injection = factor(injection, levels = c("naive_mean", "corrected_mean"))) %>%
  ggplot(., aes(analysis_type, mean, fill = injection)) +
  geom_col(position = "dodge") +
  facet_wrap(~baseline) +
  scale_fill_viridis_d(
    name = "Error assumption",
    option = "viridis",
    begin = 0,
    end = 0.9) +
  labs(
    x = "Analysis type", 
    y = "Absolute mean",
    title = "Absolute mean values at each baseline, per analysis") +
  theme_bw()+
  theme(legend.position = "top") 

all_sensitivities %>%
  dplyr::select(analysis_type, naive_mean, corrected_mean, baseline, auc) %>%
  tidyr::pivot_longer(
    cols = c(naive_mean, corrected_mean),  
    names_to = "injection",
    values_to = "mean"
  ) %>%
  dplyr::mutate(assumption = factor(case_when(
    injection == "naive_mean" ~ "Naive",
    injection == "corrected_mean" ~ "Corrected"
  ), levels = c("Naive", "Corrected"))) %>%
  dplyr::mutate(analysis_type = factor(analysis_type, levels = c("RSA","RBA_OL","RBA_CL")),
         injection = factor(injection, levels = c("naive_mean", "corrected_mean"))) %>%
  ggplot(., aes(analysis_type, mean, fill = assumption)) +
  geom_col(position = "dodge") +
  facet_grid(baseline~auc) +
  scale_fill_viridis_d(
    name = "Error assumption",
    option = "viridis",
    begin = 0,
    end = 0.9) +
  labs(
    x = "Analysis type", 
    y = "Absolute difference (delta)") +
  theme_bw()+
  theme(legend.position = "top") 

## END ##