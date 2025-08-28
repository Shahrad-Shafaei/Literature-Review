# --- Lin-Wei-Yang-Ying (LWYY) Model ---
# The LWYY model is a robust "marginal" model for total events.
# It is implemented in R using the `coxph` function with a specific
# variance estimation method. It uses the same "counting process" data
# format as the Andersen-Gill model.

# Required library
library(survival)

# --- 1. Use the Same Recurrent Event Data as the AG Model ---
# We will reuse the `recurrent_data` created in the AG model example.
# This data has start-stop intervals for each event.

set.seed(456) # Use the same seed to ensure data is identical
ids <- 1:50
recurrent_data <- data.frame()
for (i in ids) {
  n_events <- rpois(1, lambda = 1.5)
  treat_group <- sample(c("A", "B"), 1)
  if (n_events == 0) {
    recurrent_data <- rbind(recurrent_data, data.frame(
      id = i, t_start = 0, t_stop = 365, status = 0, treat = treat_group
    ))
  } else {
    event_times <- sort(runif(n_events, 1, 365))
    for (j in 1:n_events) {
      start_time <- if (j == 1) 0 else event_times[j-1]
      recurrent_data <- rbind(recurrent_data, data.frame(
        id = i, t_start = start_time, t_stop = event_times[j], status = 1, treat = treat_group
      ))
    }
    recurrent_data <- rbind(recurrent_data, data.frame(
      id = i, t_start = event_times[n_events], t_stop = 365, status = 0, treat = treat_group
    ))
  }
}

print("Using the same recurrent data as the AG model.")
head(recurrent_data)

# --- 2. Fit the LWYY Model ---
# The key difference from the AG model is the `robust` variance estimator.
# While the AG model uses `cluster(id)`, the LWYY model is specified by
# using `id = id` inside the `coxph` call to identify the patient clusters
# and setting `robust = TRUE`. This changes how the standard errors are calculated,
# making the model more robust to misspecification of the underlying event process.

# Note: The syntax can be tricky. The `id` argument links observations from the same subject.
lwyy_model <- coxph(Surv(t_start, t_stop, status) ~ treat,
                    data = recurrent_data,
                    id = id,
                    robust = TRUE)


# --- 3. View the Model Summary ---
# The output looks very similar to the AG model, but the standard errors
# (and therefore confidence intervals and p-values) will be slightly different
# due to the robust variance calculation. The HR is interpreted in the same way.

print("--- Lin-Wei-Yang-Ying (LWYY) Model Summary ---")
summary(lwyy_model)

# Compare the 'robust se' from this summary to the 'robust se' from the AG model.
# They will be different. The LWYY is generally considered a more conservative and
# reliable estimate when there is high heterogeneity.
