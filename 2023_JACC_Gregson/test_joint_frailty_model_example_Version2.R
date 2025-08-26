# Example: Joint Frailty Model for Recurrent and Terminal Events in R

# The joint frailty model assesses the effect of treatment on both types of events, accounting for shared frailty (random effect) at the subject level.
# The output will show hazard ratios for recurrent and terminal events and the frailty (association) parameter.


# Required package:
# Joint frailty models are most commonly implemented with the `frailtypack` package.
# Install if needed: install.packages("frailtypack")



library(frailtypack)

# Simulate example data for 30 subjects, each with up to 3 recurrent events and a possible terminal event (death)
set.seed(2025)
n_subjects <- 30
max_events <- 3

# Simulate a binary treatment
trt <- rbinom(n_subjects, 1, 0.5)
dat_list <- list()
death_data <- data.frame()

for (i in 1:n_subjects) {
  # Simulate recurrent event times as gap times, then cumulative
  n_ev <- sample(1:max_events, 1)
  gap_times <- rexp(n_ev, rate = 0.25 + 0.2 * trt[i])
  recur_times <- cumsum(gap_times)
  # Censor recurrent events at a random time (censoring time)
  censor_time <- runif(1, 2, 8)
  observed <- recur_times[recur_times <= censor_time]
  n_obs <- length(observed)
  if (n_obs == 0) next
  # Build data for recurrent events for this subject
  dat_list[[i]] <- data.frame(
    id = i,
    time = observed,
    event = rep(1, n_obs),
    trt = trt[i]
  )
  if (max(observed) < censor_time) {
    # Add censored recurrent event
    dat_list[[i]] <- rbind(dat_list[[i]], data.frame(
      id = i,
      time = censor_time,
      event = 0,
      trt = trt[i]
    ))
  }
  # Simulate death (terminal event) after last observed time, possibly censored
  death_time <- censor_time + rexp(1, rate = 0.10 + 0.1 * trt[i])
  death_status <- rbinom(1, 1, 0.7) # 70% chance of death observed
  death_data <- rbind(death_data, data.frame(
    id = i,
    time = death_time,
    death = death_status,
    trt = trt[i]
  ))
}

# Combine recurrent events data
recur_data <- do.call(rbind, dat_list)
# Only keep the last row for each subject for death data
death_data <- death_data[!duplicated(death_data$id), ]

# Print first few rows
cat("Recurrent events data (head):\n")
print(head(recur_data))
cat("\nTerminal event (death) data (head):\n")
print(head(death_data))

# Fit the joint frailty model
# The function jointModel() in frailtypack requires time and event for recurrent,
# and time and status for terminal (death) event, plus a subject id
jointfit <- jointModel(
  formula = Surv(time, event) ~ trt + cluster(id),
  formula.terminalEvent = Surv(time, death) ~ trt,
  data = recur_data,
  data.terminalEvent = death_data,
  n.knots = 7,     # number of knots for baseline hazard spline
  kappa = c(1, 1)  # smoothing parameters
)

# Show summary
summary(jointfit)