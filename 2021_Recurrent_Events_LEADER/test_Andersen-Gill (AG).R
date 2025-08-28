
rm(list = ls())

# Load necessary libraries
library(survival)

# Create sample data for Andersen-Gill (AG) model
# Simulate data for 3 subjects, each with multiple events (start, stop, event)

set.seed(123)
n_subjects <- 3
n_events <- 5

# Expand data so each subject can have multiple intervals
sample_data <- data.frame(
  id = rep(1:n_subjects, each = n_events),
  start = rep(seq(0, 40, by = 10)[1:n_events], n_subjects),
  stop = rep(seq(10, 50, by = 10)[1:n_events], n_subjects)
)

# Simulate event indicator (1=event occurred, 0=censored)
sample_data$event <- rbinom(n_subjects * n_events, 1, 0.6)

# Simulate a covariate (e.g., treatment group)
sample_data$trt <- rep(sample(0:1, n_subjects, replace = TRUE), each = n_events)

# Print sample data
print(sample_data)

# Fit the Andersen-Gill (AG) model using coxph with counting process style input
ag_model <- coxph(Surv(start, stop, event) ~ trt + cluster(id), data = sample_data)

# Show model summary
summary(ag_model)
