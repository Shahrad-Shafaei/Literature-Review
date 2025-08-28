# Load required package
library(survival)

# Simulate data for 5 subjects, each with up to 3 events
set.seed(2025)
n_subjects <- 5
n_events   <- 3

# Build the base data frame
dat <- data.frame(
  id = rep(1:n_subjects, each = n_events),
  event_number = rep(1:n_events, times = n_subjects)
)

# Simulate event gap times and cumulative event times for each subject
dat$gap_time <- rexp(n_subjects * n_events, rate = 0.15)
dat$stop     <- ave(dat$gap_time, dat$id, FUN = cumsum)
dat$start    <- ave(dat$stop, dat$id, FUN = function(x) c(0, head(x, -1)))

# Simulate a binary treatment covariate
dat$trt <- rep(sample(0:1, n_subjects, replace = TRUE), each = n_events)

# Event indicator: last event for each subject censored
dat$event <- 1
dat$event[seq(n_events, n_subjects*n_events, by = n_events)] <- 0

print(dat)

### 1. Wei-Lin-Weissfeld (WLW) model ###
# Each event treated as separate process, time from baseline
fit_wlw <- coxph(
  Surv(stop, event) ~ trt + strata(event_number) + cluster(id),
  data = dat
)
cat("\nWei-Lin-Weissfeld (WLW) Model Summary:\n")
print(summary(fit_wlw))

### 2. Prentice-Williams-Peterson (PWP) model ###
# Each event treated conditionally; only subjects with prior event are at risk for next
# Time is gap time since last event (start, stop)
fit_pwp <- coxph(
  Surv(start, stop, event) ~ trt + strata(event_number) + cluster(id),
  data = dat
)
cat("\nPrentice-Williams-Peterson (PWP) Model Summary:\n")
print(summary(fit_pwp))
