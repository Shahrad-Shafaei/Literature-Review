# --- Andersen-Gill (AG) Model for Total Events ---
# This model is a counting process formulation of the Cox model,
# designed to handle multiple events per subject.

# Required library
library(survival)

# --- 1. Create Sample Data for Recurrent Events ---
# This data format is different. Each row represents an interval of time for a patient.
# 'id': Patient identifier
# 't_start': The start time of the interval (e.g., 0 or time of previous event)
# 't_stop': The end time of the interval (time of event or censoring)
# 'status': 1 if an event occurred at t_stop, 0 otherwise
# 'treat': Treatment group

set.seed(456)
ids <- 1:50 # 50 patients
recurrent_data <- data.frame()
for (i in ids) {
  n_events <- rpois(1, lambda = 1.5) # Each patient has a random number of events
  treat_group <- sample(c("A", "B"), 1)
  
  if (n_events == 0) {
    # Patient with no events
    recurrent_data <- rbind(recurrent_data, data.frame(
      id = i, t_start = 0, t_stop = 365, status = 0, treat = treat_group
    ))
  } else {
    event_times <- sort(runif(n_events, 1, 365))
    # Add intervals for each event
    for (j in 1:n_events) {
      start_time <- if (j == 1) 0 else event_times[j-1]
      recurrent_data <- rbind(recurrent_data, data.frame(
        id = i, t_start = start_time, t_stop = event_times[j], status = 1, treat = treat_group
      ))
    }
    # Add final interval after the last event until end of follow-up
    recurrent_data <- rbind(recurrent_data, data.frame(
      id = i, t_start = event_times[n_events], t_stop = 365, status = 0, treat = treat_group
    ))
  }
}

# View the data for a patient with multiple events
print("Sample Recurrent Data (for patient 1):")
print(head(recurrent_data[recurrent_data$id == 1,]))


# --- 2. Fit the Andersen-Gill (AG) Model ---
# The formula is similar to Cox, but the Surv object now takes start and stop times.
# `cluster(id)` is crucial: it tells the model that multiple rows belong to the
# same patient, adjusting the standard errors to account for this correlation.
# This is what makes it a proper recurrent event model.

ag_model <- coxph(Surv(t_start, t_stop, status) ~ treat + cluster(id), data = recurrent_data)

# --- 3. View the Model Summary ---
# The interpretation is similar to the Cox model. The HR (exp(coef)) represents
# the effect of the treatment on the hazard of an event at any given time,
# considering all events.

print("--- Andersen-Gill (AG) Model Summary ---")
summary(ag_model)

# The 'robust se' in the summary indicates that the model has correctly
# accounted for the clustered nature of the data (multiple events per patient).
