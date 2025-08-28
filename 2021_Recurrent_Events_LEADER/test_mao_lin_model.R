# Load required package
library(survival)

# Simulate sample data for Mao-Lin (ML) marginal model for recurrent events
# 6 subjects, each with up to 3 events (marginal Cox model with robust SE)
set.seed(2025)
n_subjects <- 6
n_events <- 3

# Construct base data frame
dat <- data.frame(
  id = rep(1:n_subjects, each = n_events),
  event_number = rep(1:n_events, times = n_subjects)
)

# Simulate start and stop times for recurring events
dat$start <- 0
dat$gap_time <- rexp(n_subjects * n_events, rate = 0.12)
dat$stop <- ave(dat$gap_time, dat$id, FUN = cumsum)

# Binary treatment assignment per subject
dat$trt <- rep(sample(0:1, n_subjects, replace = TRUE), each = n_events)

# Event indicator (last event for each subject censored)
dat$event <- 1
dat$event[seq(n_events, n_subjects * n_events, by = n_events)] <- 0

# Print the simulated data
print(dat)

# Fit the Mao-Lin marginal model: Cox PH with robust variance for recurrent events
fit_ml <- coxph(Surv(start, stop, event) ~ trt + cluster(id), data = dat)
summary(fit_ml)

