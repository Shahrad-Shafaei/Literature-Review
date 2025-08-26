# Load required packages
library(MASS)      # for glm.nb()
library(dplyr)     # for data manipulation

# Each subject has a treatment group, follow-up time, and a simulated number of recurrent events.
# The model glm.nb(events ~ trt + offset(log(followup))) fits a Negative Binomial regression:
# events is the total event count per subject.
# trt is the treatment indicator.
# offset(log(followup)) adjusts for different follow-up times (so the model estimates rates, not just raw counts).
# exp(coef) in the output gives the incidence rate ratio (IRR) for the treatment effect.
# Interpretation: If the IRR > 1, treatment increases the event rate; if IRR < 1, it decreases the event rate.

# Simulate recurrent event data for 50 subjects
set.seed(2025)
n_subjects <- 50

# Simulate a treatment group (0/1)
trt <- rbinom(n_subjects, 1, 0.5)

# Simulate follow-up time (e.g., between 1 and 5 years)
followup <- runif(n_subjects, min = 1, max = 5)

# Simulate event rate based on treatment group
baseline_rate <- 1.2
rate <- baseline_rate * exp(0.7 * trt)   # Suppose treatment increases rate

# Simulate number of events per subject using Negative Binomial distribution
# theta (size) controls overdispersion; smaller = more overdispersion
theta <- 1.5
events <- rnbinom(n_subjects, mu = rate * followup, size = theta)

# Create data frame
dat <- data.frame(
  id = 1:n_subjects,
  trt = trt,
  followup = followup,
  events = events
)

# View the first few rows
print(head(dat))

# Fit a Negative Binomial regression model for recurrent counts
fit_nb <- glm.nb(events ~ trt + offset(log(followup)), data = dat)
summary(fit_nb)

# Interpret:
cat("\nexp(coef) gives the estimated rate ratio (incidence rate ratio) for treatment effect.\n")