# Load required package
library(survival)

# Model: The LWYY model is fit as a marginal Cox model with robust (sandwich) standard errors using cluster(id) to account for within-subject correlation.
# Interpretation: The output will show the effect of the treatment (trt) on the hazard of recurrent events, with robust inference.

# Simulate sample data for Lin-Wei-Yang-Ying (LWYY) model (marginal Cox with robust SE)
# 6 subjects, each with up to 3 recurrent events
set.seed(2025)
n_subjects <- 6
n_events <- 3

# Construct the data frame
dat <- data.frame(
  id = rep(1:n_subjects, each = n_events),
  event_number = rep(1:n_events, times = n_subjects)
)

# Simulate gap times and calculate cumulative event times per subject (for recurrent events)
dat$gap_time <- rexp(n_subjects * n_events, rate = 0.2)
dat$stop <- ave(dat$gap_time, dat$id, FUN = cumsum)
dat$start <- ave(dat$stop, dat$id, FUN = function(x) c(0, head(x, -1)))

# Simulate a binary treatment variable
dat$trt <- rep(sample(0:1, n_subjects, replace = TRUE), each = n_events)

# Simulate event indicator (last event censored for each subject)
dat$event <- 1
dat$event[seq(n_events, n_subjects * n_events, by = n_events)] <- 0

# Print the simulated data
print(dat)

# Fit the LWYY marginal model (Cox PH with robust/sandwich variance for repeated events)
fit_lwyy <- coxph(Surv(start, stop, event) ~ trt + cluster(id), data = dat)
summary(fit_lwyy)