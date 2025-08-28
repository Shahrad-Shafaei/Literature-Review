
# ===============================================================
# Recurrent Events Demo: GL and ML Curves on Toy Data (R)
# ===============================================================
# What you get:
# 1) Ghosh–Lin (GL) death-adjusted marginal mean curves via mets::recurrentMarginal
# 2) Mao–Lin (ML) proportional means for a weighted composite (recurrent + death) via Wcompo::CompoML
# No math, just plots and readable summaries.
# ---------------------------------------------------------------
# Install once (uncomment if needed)
# install.packages(c("mets", "survival", "Wcompo"))

library(mets)       # GL estimator + toy data helpers
library(survival)   # Surv(), etc.
library(Wcompo)     # ML (weighted composite) model

set.seed(20250818)

# ---------------------------------------------------------------
# 0) Simulate a toy two-arm RCT with recurrent hospitalizations + death
#    We'll make placebo a bit worse (higher recurrence and higher death).
#    Using mets::simGLcox to simulate a process consistent with GL structure.
# ---------------------------------------------------------------

n_per_arm <- 400
n <- 2 * n_per_arm
arm <- rep(c(0,1), each = n_per_arm)  # 0=Placebo, 1=Treatment

# Baseline rate functions (piecewise constant hazards) kept simple:
# rate for recurrent events > placebo arm; death hazard > placebo
# simGLcox arguments: we simulate recurrent (cause=1) and death (cause=2)
simdat <- simGLcox(n=n,
                   beta=log(0.75),     # treatment effect on recurrent events (HR ~ 0.75)
                   betad=log(0.80),    # treatment effect on death (HR ~ 0.80)
                   Z=arm,
                   max.time=36,        # 3 years of follow-up
                   rate=0.08,          # baseline recurrent hazard
                   rated=0.03)         # baseline death hazard

# simGLcox returns a long data frame with (entry, time, status, id, Z)
# status: 0=censor, 1=recurrent, 2=death
head(simdat)

# ---------------------------------------------------------------
# 1) GL curves: death-adjusted marginal mean of recurrent events
#    (expected # events per person over time, respecting that events stop after death)
# ---------------------------------------------------------------

# Two-arm estimator using strata(treatment). We map Z to factor labels.
simdat$treatment <- factor(simdat$Z, levels=c(0,1), labels=c("Placebo","Treatment"))

fit_gl <- recurrentMarginal(
  Event(entry, time, status) ~ strata(treatment) + cluster(id),
  data = simdat, cause = 1, death.code = 2
)

# Plot GL curves with simple labels
plot(fit_gl, se = FALSE, lwd = 3, col = c("black","darkgray"),
     xlab="Time (months)", ylab="Expected hospitalizations per patient",
     main = "Ghosh–Lin (GL) death-adjusted mean event curves")
legend("topleft", lwd=3, col=c("black","darkgray"), legend=levels(simdat$treatment), bty="n")

# Summaries at selected times
print(summary(fit_gl, times = c(6,12,24,36)))

# ---------------------------------------------------------------
# 2) ML curves: proportional-means model for a weighted composite
#    (recurrent events + death), using Wcompo::CompoML
#    status code for Wcompo: 1=death, 2..K=recurrent types, 0=censor
# ---------------------------------------------------------------

# Wcompo expects long-format with (id, time, status) and covariates Z.
# Our simdat already matches this; just ensure status coding matches:
# Here: status 2=death in mets; Wcompo needs death=1.
sim_ml <- simdat
sim_ml$status_ml <- sim_ml$status
sim_ml$status_ml[sim_ml$status==2] <- 1     # death -> 1
sim_ml$status_ml[sim_ml$status==1] <- 2     # recurrent -> 2

# Design matrix Z: intercept handled internally; we pass treatment as a column
Z <- model.matrix(~ treatment, data = sim_ml)[, -1, drop=FALSE]  # treatment (Treatment vs Placebo)

# Choose weights w = (w_D, w_recur). Many CV trials up-weight death; here 2:1
w <- c(2,1)

fit_ml <- CompoML(id = sim_ml$id,
                  time = sim_ml$time,
                  status = sim_ml$status_ml,
                  Z = Z,
                  w = w)

cat("\nMao–Lin (ML) proportional means for weighted composite (Death:Recur = 2:1)\n")
print(fit_ml)

# Plot model-based mean composite curves by arm
op <- par(no.readonly = TRUE)
par(mfrow=c(1,1))
plot(fit_ml, z = c(0), lwd=3, xlab="Time (months)",
     ylab="Expected weighted composite per patient",
     main="Mao–Lin (ML) weighted composite mean curves")
plot(fit_ml, z = c(1), add=TRUE, lwd=3, lty=2)
legend("topleft", lwd=3, lty=c(1,2), legend=c("Placebo","Treatment"), bty="n")
par(op)

# ---------------------------------------------------------------
# 3) Bonus: unadjusted mean curves (NA-style) for teaching contrast
#    Using reda::mcf would be another route; here we show GL vs ML only.
# ---------------------------------------------------------------

cat("\nDone. You now have:\n",
    "- GL curves (death-adjusted mean recurrent events) via mets::recurrentMarginal\n",
    "- ML curves (weighted composite mean) via Wcompo::CompoML\n")

