
rm(list = ls())


# Load required package
library(survival)

# Simulate sample data for Frailty Model
# 5 groups (clusters), each with 8 subjects
set.seed(100)
n_groups <- 5
n_per_group <- 8

# Cluster (group) variable
group <- rep(1:n_groups, each = n_per_group)

# Simulate a binary covariate (treatment)
trt <- rep(sample(0:1, n_groups, replace = TRUE), each = n_per_group)

# Simulate survival times (exponential) with random effect per group (frailty)
frailty_effect <- rep(rgamma(n_groups, shape=1, scale=1), each = n_per_group)

# let's try a frailty term (Shahrad suggestion)
# this could be smoking 
# or even gene expression level
frailty_effect <- c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5,
                    6.8, 6.8, 6.8, 6.8, 6.8, 6.8, 6.8, 6.8, 
                    1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 
                    12, 12, 12, 12, 12, 12, 12, 12, 
                    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
baseline_hazard <- 0.05
linpred <- 0.5 * trt + log(frailty_effect)
time <- rexp(n_groups * n_per_group, rate = baseline_hazard * exp(linpred))

# Right-censoring
censoring_time <- rexp(n_groups * n_per_group, rate = 0.01)
status <- as.numeric(time <= censoring_time)
obs_time <- pmin(time, censoring_time)

# Combine into a data frame
dat <- data.frame(
  id = 1:(n_groups * n_per_group),
  group = factor(group),
  trt = trt,
  time = obs_time,
  status = status
)
print(dat)

# Fit a frailty model (random effect for group)
fit_frailty <- coxph(Surv(time, status) ~ trt + frailty(group), data = dat)
summary(fit_frailty)
