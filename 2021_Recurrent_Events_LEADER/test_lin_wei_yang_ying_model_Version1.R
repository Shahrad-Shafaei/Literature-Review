# Load required package
library(survival)

# Simulate sample data for Lin-Wei-Yang-Ying (LWYY) marginal model
# 5 subjects, each with up to 3 events (marginal approach)
set.seed(2025)
n_subjects <- 5
n_events <- 3

# Build data frame
dat <- data.frame(
  id = rep(1:n_subjects, each = n_events),
  event_number = rep(1:n_events, times = n_subjects)
)

# Simulate start/stop times for recurrent events per subject
dat$start <- 0
dat$stop <- rexp(n_subjects * n_events, rate = 0.2)
dat$trt <- rep(sample(0:1, n_subjects, replace = TRUE), each = n_events)
dat$event <- 1
dat$event[seq(n_events, n_subjects * n_events, by = n_events)] <- 0  # censor last event

# Print simulated data
print(dat)

# Fit the LWYY marginal model using robust standard error
fit_lwyy <- coxph(Surv(start, stop, event) ~ trt + cluster(id), data = dat)
summary(fit_lwyy)