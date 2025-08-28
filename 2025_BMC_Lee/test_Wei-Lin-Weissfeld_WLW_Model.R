# --- Wei-Lin-Weissfeld (WLW) Model ---
# The WLW model is a "marginal" model that fits a separate Cox model for each
# event number (1st event, 2nd event, etc.) and then combines the results.
# It does not assume an order and is robust but can be less powerful.

# Required library
library(survival)

# --- 1. Create Sample Data for WLW Model ---
# The data structure is different. We need one row for each potential event
# for each patient, measuring time from the start of the study.
# 'id': Patient identifier
# 'time': Time from randomization to the k-th event
# 'status': 1 if the k-th event occurred, 0 otherwise
# 'event_num': The order of the event (1, 2, 3...)
# 'treat': Treatment group

set.seed(202)
ids <- 1:100 # 100 patients
max_events <- 3 # We will model up to the 3rd event
wlw_data <- data.frame()

for (i in ids) {
  n_events <- sample(0:max_events, 1, prob = c(0.4, 0.3, 0.2, 0.1))
  treat_group <- sample(c("A", "B"), 1)
  event_times <- sort(runif(n_events, 1, 365))
  
  for (k in 1:max_events) {
    # For each potential event number k...
    if (k <= n_events) {
      # The patient had this event
      time_k <- event_times[k]
      status_k <- 1
    } else {
      # The patient did not have this event
      # Time is the last known event time or end of follow-up
      time_k <- if(n_events > 0) event_times[n_events] else 365
      status_k <- 0
    }
    wlw_data <- rbind(wlw_data, data.frame(
      id = i, time = time_k, status = status_k, event_num = k, treat = treat_group
    ))
  }
}

# View the data for a single patient
print("Sample WLW Data (for patient with ID 1):")
print(wlw_data[wlw_data$id == 1,])


# --- 2. Fit the WLW Model ---
# The model is fit very similarly to the PWP model.
# `strata(event_num)` fits a separate model for each event number.
# `cluster(id)` combines the results and computes a robust standard error
# that accounts for the fact that the three rows for patient 1 are correlated.

wlw_model <- coxph(Surv(time, status) ~ treat + strata(event_num) + cluster(id),
                   data = wlw_data)


# --- 3. View the Model Summary ---
# The summary provides a single, overall Hazard Ratio for the treatment.
# This HR is a weighted average of the event-specific HRs.
# This approach is robust because it makes fewer assumptions about the
# relationship between events compared to the PWP model.

print("--- WLW Model Summary ---")
summary(wlw_model)
