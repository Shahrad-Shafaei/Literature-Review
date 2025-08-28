
rm(list = ls())


# Load required package
library(survival)

# Simulate sample data for Wei-Lin-Weissfeld (WLW) marginal model
# 4 subjects, each with up to 3 possible events (event order = stratum)
set.seed(2025)
n_subjects <- 4
n_events <- 3

# Create base structure
dat <- data.frame(
  id = rep(1:n_subjects, each = n_events),
  event_number = rep(1:n_events, times = n_subjects)
)

# Simulate event times (from baseline for each event, as required by WLW)
dat$time <- rexp(n_subjects * n_events, rate = 0.1)

# Simulate event indicator (last event for each subject might be censored)
dat$event <- 1
dat$event[seq(n_events, n_subjects * n_events, by = n_events)] <- 0  # Censor last event

# Simulate a binary covariate (e.g., treatment group)
dat$trt <- rep(sample(0:1, n_subjects, replace = TRUE), each = n_events)

# Print the simulated data
print(dat)

# Fit the WLW marginal model: stratify by event number (event order)
fit_wlw <- coxph(Surv(time, event) ~ trt + strata(event_number) + cluster(id), data = dat)


summary(fit_wlw)