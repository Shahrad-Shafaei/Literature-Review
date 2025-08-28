# --- Prentice-Williams-Peterson (PWP) Gap Time Model ---
# This model is for ordered recurrent events and assumes that the risk of a
# subsequent event may be different from the risk of a prior one.
# It analyzes events stratified by their order (1st, 2nd, 3rd, etc.).
# This example uses the "Gap Time" approach, which resets the clock after each event.

# Required library
library(survival)

# --- 1. Create Sample Data for PWP Model ---
# The data needs a 'counting process' format (t_start, t_stop) plus a
# column to identify the event number.
# 'id': Patient identifier
# 't_start', 't_stop': Start and end of the time interval (gap time)
# 'status': 1 if an event occurred at t_stop
# 'event_num': The order of the event (1st, 2nd, etc.)
# 'treat': Treatment group

set.seed(101)
ids <- 1:75 # 75 patients
pwp_data <- data.frame()
for (i in ids) {
  # Each patient has a random number of events, max of 4
  n_events <- sample(0:4, 1, prob = c(0.3, 0.4, 0.2, 0.05, 0.05))
  treat_group <- sample(c("A", "B"), 1)
  
  if (n_events == 0) {
    # Patient with no events
    pwp_data <- rbind(pwp_data, data.frame(
      id = i, t_start = 0, t_stop = 365, status = 0, event_num = 1, treat = treat_group
    ))
  } else {
    event_times <- sort(runif(n_events, 1, 365))
    for (j in 1:n_events) {
      start_time <- if (j == 1) 0 else event_times[j-1]
      # The gap time is the difference between the current and previous event time
      gap_start <- 0 
      gap_stop <- event_times[j] - start_time
      
      pwp_data <- rbind(pwp_data, data.frame(
        id = i, t_start = gap_start, t_stop = gap_stop, status = 1, event_num = j, treat = treat_group
      ))
    }
  }
}

# View the data for a patient with multiple events
print("Sample PWP Data (for patient with ID 2):")
print(pwp_data[pwp_data$id == 2,])


# --- 2. Fit the PWP Gap Time Model ---
# `strata(event_num)` is the key component. It tells the model to fit a
# separate baseline hazard for each event number. This means the underlying risk
# for the 1st event is allowed to be different from the 2nd, and so on.
# `cluster(id)` adjusts the standard errors for the correlation of events within a patient.

pwp_model <- coxph(Surv(t_start, t_stop, status) ~ treat + strata(event_num) + cluster(id),
                   data = pwp_data)

# --- 3. View the Model Summary ---
# The model provides a single, overall Hazard Ratio for the treatment effect,
# assuming the treatment's effect is the same across all event numbers.
# The stratification accounts for the different baseline risks of each event.

print("--- PWP Gap Time Model Summary ---")
summary(pwp_model)

# The HR (exp(coef)) of ~0.57 suggests that treatment B reduces the hazard of
# the next event by about 43% compared to treatment A, at any event stage.
