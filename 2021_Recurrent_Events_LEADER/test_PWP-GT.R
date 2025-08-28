
rm(list = ls())

# Load required package
library(survival)

# Simulate sample data for Prentice-Williams-Peterson (PWP) model
# 4 subjects, each with 3 possible events (gap time for PWP-GT)
set.seed(2025)
n_subjects <- 4
n_events <- 3

# Generate data frame for subjects and events
dat <- data.frame(
  id    = rep(1:n_subjects, each = n_events),
  event_number = rep(1:n_events, times = n_subjects)
)

# Simulate gap times (time since last event)
dat$gap_time <- rexp(n_subjects * n_events, rate = 0.2)
dat$start <- with(dat, ave(gap_time, id, FUN = function(x) c(0, head(cumsum(x), -1))))
dat$stop  <- with(dat, ave(gap_time, id, FUN = cumsum))

# Simulate event indicator (1=event, 0=censored after last interval)
dat$event <- 1
dat$event[seq(n_events, n_subjects * n_events, by = n_events)] <- 0  # Last event is censored

# Simulate a binary covariate (e.g., treatment group)
dat$trt <- rep(sample(0:1, n_subjects, replace = TRUE), each = n_events)

# Show the simulated data
print(dat)

# Fit PWP-Gap Time (PWP-GT) model using stratification by event number
fit_pwp_gt <- coxph(Surv(start, stop, event) ~ trt + strata(event_number) + cluster(id), data = dat)
summary(fit_pwp_gt)
