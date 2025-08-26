# Load required packages
library(survival)
library(ggplot2)


# This code uses the reda package for mean cumulative function (MCF) and AUC estimation for recurrent events.
# The AUC here is the area under the MCF curve, which estimates the mean total number of events per subject up to a fixed time (e.g., 5).
# Treatment group differences in AUC reflect differences in expected recurrent event burden.
# Install reda with install.packages("reda") if needed.
# The plot and printed AUC values give a practical summary of recurrent event experience by group.


# Simulate recurrent event data: 10 subjects, up to 4 events each
set.seed(123)
n_subjects <- 10
n_events <- 4

dat <- data.frame(
  id = rep(1:n_subjects, each = n_events),
  event_number = rep(1:n_events, times = n_subjects)
)

# Simulate event times as cumulative sum of gap times (exponential)
dat$gap_time <- rexp(n_subjects * n_events, rate = 0.5)
dat$event_time <- ave(dat$gap_time, dat$id, FUN = cumsum)

# Simulate event indicator (last event censored)
dat$event <- 1
dat$event[seq(n_events, n_subjects * n_events, by = n_events)] <- 0

# Simulate treatment group
dat$trt <- rep(sample(0:1, n_subjects, replace = TRUE), each = n_events)

# Print sample data
print(dat)

# Prepare data for mean cumulative function (MCF) estimation
# We'll estimate the mean cumulative number of events (Nelson-Aalen estimator)
library(reda) # For AUC and MCF, install if needed: install.packages("reda")

# Prepare data in counting process format
dat_recur <- dataRecur(
  id = dat$id,
  time = dat$event_time,
  event = dat$event,
  group = dat$trt
)

# Estimate mean cumulative function (MCF) for each treatment group
mcf_fit <- mcf(Recur(event_time, id, event) ~ trt, data = dat)

# Plot MCFs
autoplot(mcf_fit) + ggtitle("Mean Cumulative Function (MCF) by Treatment Group")

# Calculate area under the MCF curve (AUC) up to a fixed time (e.g., t = 5)
max_time <- 5
mcf_est <- summary(mcf_fit, times = seq(0, max_time, by = 0.1))

# Calculate AUC (trapezoidal rule) for each group
auc_table <- by(mcf_est$MCF, mcf_est$trt, function(mcf_group) {
  time <- mcf_est$time[mcf_est$trt == as.numeric(names(mcf_group)[1])]
  auc <- sum(diff(time) * (head(mcf_group, -1) + tail(mcf_group, -1)) / 2)
  return(auc)
})

print("Estimated Area Under the MCF Curve (AUC) up to time 5 for each group:")
print(auc_table)