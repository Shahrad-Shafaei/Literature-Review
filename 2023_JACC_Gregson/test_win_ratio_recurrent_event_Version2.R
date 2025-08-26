# Install and load required packages
# install.packages("WR")
library(WR)


# Data are simulated for two groups (treatment and control), each with several subjects and multiple recurrent events.
# Each subject’s primary outcome is the number of events (fewer is better), and the secondary is time to terminal event (longer is better).
# The wr() function from the WR package calculates the Win Ratio: for all possible pairs (one from each group), it counts how often a subject in the treatment group "wins" (has fewer events, or if tied, survives longer) vs. the control group.
# The code prints both the event-level and subject-level data and the Win Ratio results.
# The code assumes Win Ratio is based on:
# Number of recurrent events (lower is better)
# If tied, longer survival/censoring time is better
# Note:
  # Win Ratio is more common for time-to-first-event or composite endpoints, but the above approach adapts it to recurrent events (using total event count as primary outcome), as seen in some recent methodological papers.
  # For more complex recurrent event Win Ratio (e.g., pairwise comparison of each event), more advanced code or packages may be required.                                                                                                                                                                                                                                                             
                                                                                                                                                                           


# Simulate sample data for Win Ratio analysis in recurrent events context
# 10 subjects per group, 2 groups (trt = 0 or 1)
set.seed(2025)
n_per_group <- 10
n_subjects <- 2 * n_per_group

# Simulate treatment assignment
trt <- rep(0:1, each = n_per_group)

# Simulate number of recurrent events per subject (higher for trt=0)
n_events <- rpois(n_subjects, lambda = ifelse(trt == 0, 3, 2))

# Simulate recurrent event times (from baseline)
event_times <- lapply(n_events, function(k) if (k > 0) sort(runif(k, 0, 10)) else numeric(0))

# Simulate terminal (censoring) time
terminal_time <- runif(n_subjects, 8, 12)

# Build a long-format data frame: one row per subject per event
df_events <- do.call(rbind, lapply(1:n_subjects, function(i) {
  if (length(event_times[[i]]) == 0) return(NULL)
  data.frame(
    id = i,
    trt = trt[i],
    event_time = event_times[[i]],
    terminal_time = terminal_time[i]
  )
}))

# Print first few rows of recurrent event data
print(head(df_events))

# Prepare data for WR::wr() function
# The WR package expects a data frame with subject-level summary
df_summary <- data.frame(
  id = 1:n_subjects,
  trt = trt,
  n_events = sapply(event_times, length),
  terminal_time = terminal_time
)

# Print subject-level summary
print(df_summary)

# Use the Win Ratio (WR) comparing treatment (trt=1) vs control (trt=0)
# Primary: number of recurrent events (lower is better)
# Secondary: time to terminal event (later is better)

# The WR package's wr() function expects:
# - group: group indicator (0/1)
# - primary: numeric vector (e.g., number of events)
# - secondary: numeric vector (e.g., terminal time)

wr_result <- wr(
  group = df_summary$trt,
  primary = -df_summary$n_events,    # Negative so that fewer events is considered a 'win'
  secondary = df_summary$terminal_time
)

print(wr_result)