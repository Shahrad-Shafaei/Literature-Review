# Load required package
library(survival)

# Simulate sample data for the Ghosh-Lin (GL) marginal model
# 5 subjects, each with up to 3 recurrent events (GL is a marginal model for recurrent events and terminal event)
set.seed(2025)
n_subjects <- 5
n_events <- 3

# Construct a data frame for recurrent events
dat <- data.frame(
  id = rep(1:n_subjects, each = n_events),
  event_number = rep(1:n_events, times = n_subjects)
)

# Simulate event times (gap times, then cumulative sum for each subject)
gap_times <- rexp(n_subjects * n_events, rate = 0.15)
dat$time <- ave(gap_times, dat$id, FUN = cumsum)
dat$start <- 0
dat$stop <- dat$time

# Simulate event indicator (last event is censored for each subject)
dat$event <- 1
dat$event[seq(n_events, n_subjects * n_events, by = n_events)] <- 0

# Simulate a binary treatment covariate
dat$trt <- rep(sample(0:1, n_subjects, replace = TRUE), each = n_events)

# Print the simulated data
print(dat)

# Fit the Ghosh-Lin marginal model using robust variance (sandwich estimator)
# This is essentially a Cox proportional hazards model with robust SE for recurrent events
fit_gl <- coxph(Surv(start, stop, event) ~ trt + cluster(id), data = dat)
summary(fit_gl)