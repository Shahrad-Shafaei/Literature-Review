
rm(list = ls())

# Load required packages
# Install if not already installed: install.packages("frailtypack")
library(frailtypack)

# Simulate sample data for Joint Frailty Model (recurrent events + terminal event)
set.seed(1234)
n_subjects <- 100

# Simulate baseline covariate
covar <- rbinom(n_subjects, 1, 0.5)

# Simulate frailty for each subject
frailty <- rgamma(n_subjects, shape = 1, scale = 1)

# Simulate number of recurrent events per subject (Poisson with frailty)
n_events <- rpois(n_subjects, lambda = 2 * frailty)

# For each recurrent event, simulate event times (gap time, exponential)
recurrent_data <- do.call(rbind, lapply(1:n_subjects, function(i) {
  if(n_events[i] == 0) return(NULL)
  gap_times <- rexp(n_events[i], rate = 0.3 * exp(0.5 * covar[i]) * frailty[i])
  times <- cumsum(gap_times)
  data.frame(
    id = i,
    time = times,
    event = 1,
    covar = covar[i]
  )
}))
# Add a terminal event (death) for each subject
death_time <- rexp(n_subjects, rate = 0.05 * exp(0.8 * covar) * frailty)
death_event <- rbinom(n_subjects, 1, 0.6) # 60% die during study

terminal_data <- data.frame(
  id = 1:n_subjects,
  time = death_time,
  event = death_event,
  covar = covar
)

# Combine recurrent and terminal events for joint frailty model
# 'recurrent_data' for recurrent event; 'terminal_data' for terminal event
# The package requires a specific format:
#   recurrent: id, time, event, covariate(s)
#   terminal: id, time, event, covariate(s)

?frailtyPenal
# Fit the joint frailty model
fit_joint <- frailtyPenal(
  Surv(time, event) ~ covar + cluster(id),
  formula.terminalEvent = ~ covar,
  recurrentData = recurrent_data,
  terminalEvent = terminal_data,
  data = recurrent_data,
  data.terminalEvent = terminal_data,
  n.knots = 7,
  kappa = c(1, 1)
)

# Print summary of the fitted joint frailty model
summary(fit_joint)
