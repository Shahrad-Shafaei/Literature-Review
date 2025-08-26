# Title: R Code to Reproduce Table 2 from Bhatt and Mehta (2016) NEJM
# Description: This script simulates an adaptive clinical trial with sample size
#              re-estimation based on an interim analysis, as described in the
#              CHAMPION PHOENIX trial example. (Corrected Version)

# Load necessary libraries
# No special libraries are needed for this simulation, we'll use base R.

# --- Simulation Parameters ---

# Set a seed for reproducibility of the simulation results
set.seed(123)

# Number of simulations to run for each scenario
n_simulations <- 100000

# Overall one-sided significance level (alpha)
alpha <- 0.025

# Initial planned maximum sample size (N_max)
n_initial <- 10900

# Timing of the interim analysis (as a fraction of n_initial)
interim_fraction <- 0.7
n_interim <- n_initial * interim_fraction # 7630

# --- Efficacy Boundaries ---

# O'Brien-Fleming-like gamma(-5) spending function approximated for two looks
# Critical value for the interim analysis (for early stopping for efficacy)
# This corresponds to the boundary for the "Favorable" zone
z_alpha_interim <- 2.797 # Corresponds to a very small alpha spent

# Critical value for the final analysis in a non-adaptive design
z_alpha_final_nonadaptive <- 1.98 # Adjusted for the interim look

# --- Adaptation Rule Parameters ---

# Define the zones based on observed relative risk reduction (RRR)
# Unfavorable zone: RRR < 13.6%
# Promising zone: 13.6% <= RRR <= 21.2%
# Favorable zone: RRR > 21.2%
promising_zone_lower_rrr <- 0.136
promising_zone_upper_rrr <- 0.212

# Target conditional power if results fall in the promising zone
target_conditional_power <- 0.90

# --- Function to run the simulation for one scenario ---

simulate_trial <- function(p_control, true_rrr) {
  p_experimental <- p_control * (1 - true_rrr)
  
  # Store results
  outcomes <- data.frame(
    zone = character(n_simulations),
    final_z_adaptive = numeric(n_simulations),
    final_z_nonadaptive = numeric(n_simulations),
    sample_size_adaptive = numeric(n_simulations),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_simulations) {
    # --- Stage 1: Interim Analysis ---
    
    # Simulate events for the first stage (interim)
    events_control_s1 <- rbinom(1, n_interim / 2, p_control)
    events_experimental_s1 <- rbinom(1, n_interim / 2, p_experimental)
    
    # Observed event rates at interim
    p_hat_control_s1 <- events_control_s1 / (n_interim / 2)
    p_hat_experimental_s1 <- events_experimental_s1 / (n_interim / 2)
    
    # Pooled probability under the null hypothesis
    p_hat_pooled_s1 <- (events_control_s1 + events_experimental_s1) / n_interim
    
    # Add defensive check for NA values before the 'if' condition
    if (is.na(p_hat_pooled_s1)) p_hat_pooled_s1 <- 0
    
    # Avoid division by zero if p_hat_pooled_s1 is 0 or 1
    if (p_hat_pooled_s1 == 0 || p_hat_pooled_s1 == 1) {
      z1 <- 0
    } else {
      # Calculate Z-statistic at interim
      z1 <- (p_hat_control_s1 - p_hat_experimental_s1) /
        sqrt(p_hat_pooled_s1 * (1 - p_hat_pooled_s1) * (4 / n_interim))
    }
    
    # Observed RRR at interim
    observed_rrr_s1 <- 1 - (p_hat_experimental_s1 / p_hat_control_s1)
    
    # FIX: Check for NA, NaN, and Inf to prevent crash in the next 'if' statement
    if(is.na(observed_rrr_s1) || is.nan(observed_rrr_s1) || is.infinite(observed_rrr_s1)) {
      observed_rrr_s1 = 0
    }
    
    # --- Adaptation Decision ---
    n_final_adaptive <- n_initial
    
    if (observed_rrr_s1 < promising_zone_lower_rrr) {
      outcomes$zone[i] <- "Unfavorable"
    } else if (observed_rrr_s1 <= promising_zone_upper_rrr) {
      outcomes$zone[i] <- "Promising"
      # Recalculate sample size for target conditional power
      delta_hat <- log(1 - observed_rrr_s1)
      
      # Handle case where delta_hat is infinite
      if (is.infinite(delta_hat) || is.na(delta_hat) || delta_hat == 0) {
        n2_required <- 0
      } else {
        # Required sample size for the second stage
        n2_required <- ( (qnorm(1-alpha) - qnorm(1-target_conditional_power))^2 ) / (delta_hat^2 / (4/p_control)) - n_interim
      }
      
      n2_required <- max(0, n2_required)
      
      n_final_adaptive <- n_interim + n2_required
      # Cap the sample size to a reasonable maximum
      n_final_adaptive <- min(n_final_adaptive, 20000)
      
    } else {
      outcomes$zone[i] <- "Favorable"
      # If z1 crosses the efficacy boundary, we could stop early
      if (z1 >= z_alpha_interim) {
        n_final_adaptive <- n_interim
      }
    }
    
    outcomes$sample_size_adaptive[i] <- n_final_adaptive
    
    # --- Stage 2: Final Analysis ---
    
    # Non-adaptive path
    n2_nonadaptive <- n_initial - n_interim
    events_control_s2_na <- rbinom(1, n2_nonadaptive / 2, p_control)
    events_exp_s2_na <- rbinom(1, n2_nonadaptive / 2, p_experimental)
    
    total_events_control_na <- events_control_s1 + events_control_s2_na
    total_events_exp_na <- events_experimental_s1 + events_exp_s2_na
    
    p_hat_control_final_na <- total_events_control_na / (n_initial / 2)
    p_hat_exp_final_na <- total_events_exp_na / (n_initial / 2)
    p_hat_pooled_final_na <- (total_events_control_na + total_events_exp_na) / n_initial
    
    if (is.na(p_hat_pooled_final_na)) p_hat_pooled_final_na <- 0 # Defensive check
    
    if(p_hat_pooled_final_na == 0 || p_hat_pooled_final_na == 1) {
      outcomes$final_z_nonadaptive[i] <- 0
    } else {
      outcomes$final_z_nonadaptive[i] <- (p_hat_control_final_na - p_hat_exp_final_na) /
        sqrt(p_hat_pooled_final_na * (1 - p_hat_pooled_final_na) * (4 / n_initial))
    }
    
    # Adaptive path
    n2_adaptive <- n_final_adaptive - n_interim
    if (n2_adaptive > 0) {
      events_control_s2_a <- rbinom(1, n2_adaptive / 2, p_control)
      events_exp_s2_a <- rbinom(1, n2_adaptive / 2, p_experimental)
      
      total_events_control_a <- events_control_s1 + events_control_s2_a
      total_events_exp_a <- events_experimental_s1 + events_exp_s2_a
    } else { # Early stopping case
      total_events_control_a <- events_control_s1
      total_events_exp_a <- events_experimental_s1
    }
    
    p_hat_control_final_a <- total_events_control_a / (n_final_adaptive / 2)
    p_hat_exp_final_a <- total_events_exp_a / (n_final_adaptive / 2)
    p_hat_pooled_final_a <- (total_events_control_a + total_events_exp_a) / n_final_adaptive
    
    if (is.na(p_hat_pooled_final_a)) p_hat_pooled_final_a <- 0 # Defensive check
    
    if(p_hat_pooled_final_a == 0 || p_hat_pooled_final_a == 1 || n_final_adaptive == 0) {
      outcomes$final_z_adaptive[i] <- 0
    } else {
      z2 <- (p_hat_control_final_a - p_hat_exp_final_a) /
        sqrt(p_hat_pooled_final_a * (1-p_hat_pooled_final_a) * (4/n_final_adaptive))
      outcomes$final_z_adaptive[i] <- z2
    }
  }
  
  return(outcomes)
}

# --- Function to analyze and display results ---

analyze_results <- function(results, p_control, true_rrr) {
  # Overall Power (Unconditional)
  power_nonadaptive <- mean(results$final_z_nonadaptive >= z_alpha_final_nonadaptive, na.rm = TRUE)
  power_adaptive <- mean(results$final_z_adaptive >= z_alpha_final_nonadaptive, na.rm = TRUE)
  
  # Average Sample Size
  avg_n_adaptive <- mean(results$sample_size_adaptive, na.rm = TRUE)
  
  # Probabilities of entering each zone
  prob_unfavorable <- mean(results$zone == "Unfavorable", na.rm = TRUE)
  prob_promising <- mean(results$zone == "Promising", na.rm = TRUE)
  prob_favorable <- mean(results$zone == "Favorable", na.rm = TRUE)
  
  # --- Robust calculation of Conditional Power and Average N ---
  # Helper function to safely calculate conditional power
  get_cp <- function(zone_name, z_scores) {
    subset <- z_scores[results$zone == zone_name]
    if (length(subset) == 0) return(0) # Return 0 if zone was never entered
    mean(subset >= z_alpha_final_nonadaptive, na.rm = TRUE)
  }
  
  # Helper function to safely calculate average sample size
  get_avg_n <- function(zone_name) {
    subset <- results$sample_size_adaptive[results$zone == zone_name]
    if (length(subset) == 0) return(NA) # Return NA if zone was never entered
    mean(subset, na.rm = TRUE)
  }
  
  # Conditional Power
  cp_nonadaptive_unfavorable <- get_cp("Unfavorable", results$final_z_nonadaptive)
  cp_adaptive_unfavorable <- get_cp("Unfavorable", results$final_z_adaptive)
  
  cp_nonadaptive_promising <- get_cp("Promising", results$final_z_nonadaptive)
  cp_adaptive_promising <- get_cp("Promising", results$final_z_adaptive)
  
  cp_nonadaptive_favorable <- get_cp("Favorable", results$final_z_nonadaptive)
  cp_adaptive_favorable <- get_cp("Favorable", results$final_z_adaptive)
  
  # Average sample size within each zone
  avg_n_unfavorable <- get_avg_n("Unfavorable")
  avg_n_promising <- get_avg_n("Promising")
  avg_n_favorable <- get_avg_n("Favorable")
  
  # Print results in a table format
  cat(sprintf("\n--- Results for Control Rate: %.3f, True RRR: %.2f ---\n", p_control, true_rrr))
  cat(sprintf("Overall Power (Non-Adaptive): %.2f\n", power_nonadaptive))
  cat(sprintf("Overall Power (Adaptive):     %.2f\n", power_adaptive))
  cat(sprintf("Avg. Sample Size (Adaptive):  %.0f\n\n", avg_n_adaptive))
  
  cat("Zone          | Prob. Enter | Cond. Power (A) | Cond. Power (NA) | Avg. N (A) | Avg. N (NA)\n")
  cat("-----------------------------------------------------------------------------------------\n")
  cat(sprintf("Unfavorable   | %.2f        | %.2f            | %.2f             | %-10.0f | %-10.0f\n",
              prob_unfavorable, cp_adaptive_unfavorable, cp_nonadaptive_unfavorable, avg_n_unfavorable, n_initial))
  cat(sprintf("Promising     | %.2f        | %.2f            | %.2f             | %-10.0f | %-10.0f\n",
              prob_promising, cp_adaptive_promising, cp_nonadaptive_promising, avg_n_promising, n_initial))
  cat(sprintf("Favorable     | %.2f        | %.2f            | %.2f             | %-10.0f | %-10.0f\n",
              prob_favorable, cp_adaptive_favorable, cp_nonadaptive_favorable, avg_n_favorable, n_initial))
}


# --- Run Scenarios from Table 2 ---

# Scenario 1: p_control = 5.1%, RRR = 24%, 21%, 18%
p_control_1 <- 0.051
rrr_scenarios_1 <- c(0.24, 0.21, 0.18)
for (rrr in rrr_scenarios_1) {
  results <- simulate_trial(p_control_1, rrr)
  analyze_results(results, p_control_1, rrr)
}

# Scenario 2: p_control = 4.75%, RRR = 24%, 21%, 18%
p_control_2 <- 0.0475
rrr_scenarios_2 <- c(0.24, 0.21, 0.18)
for (rrr in rrr_scenarios_2) {
  results <- simulate_trial(p_control_2, rrr)
  analyze_results(results, p_control_2, rrr)
}
