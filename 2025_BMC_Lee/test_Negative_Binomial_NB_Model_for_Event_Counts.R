# --- Negative Binomial (NB) Model for Event Counts ---
# This model is used to analyze event rates. It doesn't use time-to-event data directly,
# but rather the total count of events for each patient over a follow-up period.
# It's particularly good when there is high heterogeneity (overdispersion).

# Required library for NB regression
library(MASS)

# --- 1. Create Sample Data for Event Counts ---
# We need one row per patient.
# 'id': Patient identifier
# 'event_count': Total number of events for that patient
# 'follow_up_time': Total time the patient was observed
# 'treat': Treatment group

set.seed(789)
n_patients <- 150
count_data <- data.frame(
  id = 1:n_patients,
  event_count = rnbinom(n_patients, size = 0.5, mu = 3), # `size` is the dispersion parameter
  follow_up_time = rep(365, n_patients), # Assume everyone followed for 1 year
  treat = sample(c("A", "B"), n_patients, replace = TRUE)
)

# View the first few rows
print("Sample Data for Event Count Analysis:")
head(count_data)

# --- 2. Fit the Negative Binomial Model ---
# We model the `event_count` based on the treatment group.
# `offset(log(follow_up_time))` is important. It accounts for different
# follow-up times, effectively turning the model from one of counts to one of *rates*.

nb_model <- glm.nb(event_count ~ treat + offset(log(follow_up_time)), data = count_data)

# --- 3. View the Model Summary ---
# The summary gives us coefficients on the log scale.
# To get the Rate Ratio (RR), we need to exponentiate the coefficient for 'treatB'.
# - RR < 1 means treatment B is associated with a lower event rate.
# - RR > 1 means treatment B is associated with a higher event rate.

print("--- Negative Binomial (NB) Model Summary ---")
summary(nb_model)

# --- 4. Calculate the Rate Ratio (RR) ---
# Extract the coefficient for the treatment and its confidence interval
coef_treat <- summary(nb_model)$coefficients["treatB", "Estimate"]
se_treat <- summary(nb_model)$coefficients["treatB", "Std. Error"]

rate_ratio <- exp(coef_treat)
ci_lower <- exp(coef_treat - 1.96 * se_treat)
ci_upper <- exp(coef_treat + 1.96 * se_treat)

print(paste0("Rate Ratio (RR) for Treat B vs. A: ", round(rate_ratio, 2)))
print(paste0("95% CI: [", round(ci_lower, 2), ", ", round(ci_upper, 2), "]"))

# An RR of 0.74 means the event rate in group B is 26% lower than in group A.
