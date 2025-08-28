
# =====================================================================
# Recurrent Events Demo: GL & ML Curves on Toy Data + Real Dataset
# =====================================================================
# Part A: Simulated two-arm trial (from prior script)
# Part B: Real recurrent-events dataset (frailtypack::readmission)
#         We create a synthetic two-arm split (for teaching only) so
#         curves differ by arm while using *real* event/death patterns.
#
# Outputs (both parts):
#  - GL curves (death-adjusted marginal mean recurrent events)
#  - ML curves (weighted composite of death + recurrent events)
#
# Tweakable knobs:
#   WEIGHTS <- c(3, 1)  # death : recurrent weighting for ML curves
#   ARM_RULE <- "random" # how to split into two arms: "random" or "by_covariate"
#   COV_NAME <- "sex"    # if ARM_RULE == "by_covariate", use this column to split
# =====================================================================

# ---- Install (uncomment if needed) ----
# install.packages(c("mets", "survival", "Wcompo", "frailtypack"))

rm(list = ls())


library(mets)        # GL estimator
library(survival)    # Surv()
library(Wcompo)      # ML estimator
library(frailtypack) # real dataset: readmission
set.seed(20250818)

WEIGHTS  <- c(3, 1)        # Death:Recurrent weight for ML
ARM_RULE <- "random"       # "random" or "by_covariate"
COV_NAME <- "sex"          # used if ARM_RULE == "by_covariate"

# ---------------------------------------------------------------------
# Helper: pretty plotting without fancy dependencies
pp_legend <- function(pos="topleft", labs, lwd=3, lty=c(1,2), bty="n") {
  legend(pos, legend=labs, lwd=lwd, lty=lty, bty=bty)
}

# =====================================================================
# Part A) Simulated two-arm trial (reprise)
# =====================================================================
message("\n=== Part A: Simulation-based GL/ML curves (as before) ===")

n_per_arm <- 400
n <- 2 * n_per_arm
arm <- rep(c(0,1), each = n_per_arm)  # 0=Placebo, 1=Treatment

?simGLcox

simdat <- simGLcox(n=n,
                   beta=log(0.75),     # treatment effect on recurrent events
                   betad=log(0.80),    # treatment effect on death
                   Z=arm,
                   max.time=36,        # months
                   rate=0.08,
                   rated=0.03)

simdat$treatment <- factor(simdat$Z, levels=c(0,1), labels=c("Placebo","Treatment"))

# GL
fit_gl_sim <- recurrentMarginal(
  Event(entry, time, status) ~ strata(treatment) + cluster(id),
  data = simdat, cause = 1, death.code = 2
)
plot(fit_gl_sim, se=FALSE, lwd=3, col=c("black","darkgray"),
     xlab="Time (months)", ylab="Expected recurrent events per patient",
     main="GL (death-adjusted) — Simulation")
pp_legend(labs=levels(simdat$treatment))

# ML
sim_ml <- within(simdat, {
  status_ml <- status
  status_ml[status==2] <- 1  # death -> 1
  status_ml[status==1] <- 2  # recurrent -> 2
})
Z_sim <- model.matrix(~ treatment, data=sim_ml)[, -1, drop=FALSE]
fit_ml_sim <- CompoML(id=sim_ml$id, time=sim_ml$time, status=sim_ml$status_ml,
                      Z=Z_sim, w=WEIGHTS)

plot(fit_ml_sim, z=c(0), lwd=3, xlab="Time (months)",
     ylab="Expected weighted composite per patient",
     main=paste0("ML (weighted composite, w=", WEIGHTS[1], ":", WEIGHTS[2], ") — Simulation"))
plot(fit_ml_sim, z=c(1), add=TRUE, lwd=3, lty=2)
pp_legend(labs=c("Placebo","Treatment"))

cat("\n[Simulation] ML summary:\n"); print(fit_ml_sim)

# =====================================================================
# Part B) Real dataset: frailtypack::readmission
# =====================================================================
message("\n=== Part B: Real dataset (frailtypack::readmission) ===")

data(readmission)
df <- readmission
cat("\nColumns in readmission:\n"); print(names(df)); cat("\nHead:\n"); print(head(df, 3))

# The dataset provides recurrent hospitalizations and a terminal event (death).
# Common columns include:
#   'id'          : subject identifier
#   'time'        : time of recurrent event or censoring (months)
#   'event'       : 1 if recurrent hospitalization occurred, 0 otherwise
#   'time2' or 'timeDeath' : time to death or censoring for terminal event
#   'death'       : 1 if death occurred, 0 otherwise
# There may be additional covariates (e.g., sex, age, dukes). Names can vary by version.

# --- Identify columns robustly ---
find_col <- function(cands) {
  ix <- which(tolower(names(df)) %in% tolower(cands))
  if (length(ix) == 0) return(NA) else return(names(df)[ix[1]])
}

col_id    <- find_col(c("id","ID","patient","subject"))
col_timeR <- find_col(c("time","timeR","time_recur","tstop","stop"))
col_evtR  <- find_col(c("event","recur","statusR","eventR"))
col_timeD <- find_col(c("time2","timeDeath","tdeath","time_death","tdeath2"))
col_evtD  <- find_col(c("death","statusD","eventD","D","terminal"))

mandatory <- c(col_id, col_timeR, col_evtR, col_timeD, col_evtD)
if (any(is.na(mandatory))) {
  stop("Could not locate required columns. Please open the dataset and set col_* manually at the top of Part B.")
}

# --- Build a long-format event history with both recurrent and death ---
# We'll create one row per recurrent event and one (optional) row for death.
# For GL: status codes 1=recurrent, 2=death
# For ML: 1=death, 2=recurrent

# Keep only columns we need
x <- df[, c(col_id, col_timeR, col_evtR, col_timeD, col_evtD)]
names(x) <- c("id","timeR","eventR","timeD","death")

# Recurrent rows
rec <- subset(x, eventR == 1, select=c("id","timeR"))
rec$status <- 1  # recurrent
names(rec)[2] <- "time"

# Death rows (one per subject with death==1)
dea <- unique(x[, c("id","timeD","death")])
dea <- subset(dea, death == 1)
names(dea)[2] <- "time"
dea$status <- 2  # death
dea$death <- NULL

# Combine
long <- rbind(
  rec[, c("id","time","status")],
  dea[, c("id","time","status")]
)
# Remove missing/zero/negative times if any
long <- long[is.finite(long$time) & long$time > 0, ]
long <- long[order(long$id, long$time, long$status), ]

# ---- Create a teaching two-arm split on *real* event history ----
# Option A: Random 1:1 split
if (ARM_RULE == "random") {
  u <- aggregate(time ~ id, data=long, FUN=length)  # one row per id
  set.seed(20250818)
  assign <- sample(c(0,1), size=nrow(u), replace=TRUE)
  arm_map <- setNames(assign, u$id)
  long$treatment <- factor(ifelse(arm_map[as.character(long$id)]==1, "Treatment", "Placebo"))
} else {
  # Option B: by a covariate threshold (e.g., sex, age)
  if (!(COV_NAME %in% names(df))) stop("COV_NAME not found in dataset.")
  subj <- unique(df[, c(col_id, COV_NAME)])
  names(subj) <- c("id","cov")
  # Split by median for numeric, by first level for factor
  if (is.numeric(subj$cov)) {
    cutp <- stats::median(subj$cov, na.rm=TRUE)
    subj$arm <- ifelse(subj$cov > cutp, "Treatment", "Placebo")
  } else {
    lev <- levels(factor(subj$cov))
    subj$arm <- ifelse(subj$cov == lev[1], "Placebo", "Treatment")
  }
  long <- merge(long, subj[, c("id","arm")], by="id", all.x=TRUE)
  long$treatment <- factor(long$arm, levels=c("Placebo","Treatment"))
  long$arm <- NULL
}

# ---- GL curves on real data ----
# Need entry times; we use counting-process style with entry=0 (calendar-time version).
# For pedagogy this is fine; for analysis you may prefer start-stop intervals.
long$entry <- 0
fit_gl_real <- recurrentMarginal(
  Event(entry, time, status) ~ strata(treatment) + cluster(id),
  data = long, cause = 1, death.code = 2
)

plot(fit_gl_real, se=FALSE, lwd=3, col=c("black","darkgray"),
     xlab="Time", ylab="Expected recurrent events per patient",
     main="GL (death-adjusted) — Real dataset")
pp_legend(labs=levels(long$treatment))

cat("\n[Real] GL summary at selected times:\n")
print(summary(fit_gl_real, times = quantile(long$time, probs=c(0.25,0.5,0.75), na.rm=TRUE)))

# ---- ML curves on real data ----
real_ml <- long
real_ml$status_ml <- ifelse(real_ml$status==2, 1, 2)  # death->1, recurrent->2
Z_real <- model.matrix(~ treatment, data=real_ml)[, -1, drop=FALSE]

fit_ml_real <- CompoML(id=real_ml$id, time=real_ml$time, status=real_ml$status_ml,
                       Z=Z_real, w=WEIGHTS)

plot(fit_ml_real, z=c(0), lwd=3, xlab="Time",
     ylab=paste0("Expected weighted composite per patient (w=", WEIGHTS[1], ":", WEIGHTS[2], ")"),
     main="ML (weighted composite) — Real dataset")
plot(fit_ml_real, z=c(1), add=TRUE, lwd=3, lty=2)
pp_legend(labs=c("Placebo","Treatment"))

cat("\n[Real] ML summary:\n"); print(fit_ml_real)

cat("\nDone.\n- GL curves reflect death-adjusted mean recurrent events.\n",
    "- ML curves reflect a weighted composite of death and recurrence (weights set at ",
    WEIGHTS[1], ":", WEIGHTS[2], ").\n", sep="")
