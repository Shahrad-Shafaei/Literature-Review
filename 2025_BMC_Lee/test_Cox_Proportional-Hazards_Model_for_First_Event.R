# --- Cox Proportional-Hazards Model for First Event ---
# This model is the standard for "time-to-first-event" analysis.
# It only considers the first event for each patient.

# Required library for survival analysis
library(survival)

# --- 1. Create Sample Data ---
# We'll create a dataset with 100 patients.
# 'id': Patient identifier
# 'time': Time (in days) to an event or end of follow-up
# 'status': 1 if an event occurred, 0 if not (censored)
# 'treat': Treatment group ('A' or 'B')

set.seed(123) # for reproducibility
n_patients <- 100
first_event_data <- data.frame(
  id = 1:n_patients,
  time = round(runif(n_patients, 50, 500)),
  status = rbinom(n_patients, 1, 0.4), # ~40% of patients have an event
  treat = sample(c("A", "B"), n_patients, replace = TRUE)
)

# View the first few rows of the data
print("Sample Data for First-Event Analysis:")
head(first_event_data)

# --- 2. Fit the Cox Model ---
# We are modeling the time to an event (`time`, `status`) based on the treatment group (`treat`).
# `Surv(time, status)` creates the survival object.
# `~ treat` specifies that treatment is the predictor variable.

cox_model <- coxph(Surv(time, status) ~ treat, data = first_event_data)

# --- 3. View the Model Summary ---
# The summary provides the Hazard Ratio (HR), confidence intervals, and p-value.
# - 'coef': The log hazard ratio. A negative value suggests the treatment is protective.
# - 'exp(coef)': This is the Hazard Ratio (HR).
#   - HR < 1 suggests treatment 'B' reduces the hazard (risk) of the event compared to 'A'.
#   - HR > 1 suggests treatment 'B' increases the hazard.
# - 'Pr(>|z|)': The p-value. If < 0.05, the effect is statistically significant.

print("--- Cox Model Summary ---")
summary(cox_model)

# The output shows an HR (exp(coef)) of ~0.69. This suggests that
# the hazard of an event in group B is about 31% lower than in group A.
# The p-value (0.16) is not significant, so we cannot conclude this is a true effect.
