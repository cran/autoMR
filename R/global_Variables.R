# ==============================================================================
# Global variables declaration for R CMD check
# Centralized bindings for:
#   (1) dplyr/ggplot2 NSE column names
#   (2) JAGS model symbols referenced in model functions/strings
# ==============================================================================

utils::globalVariables(c(
  # ---- Common MR input/output columns (NSE in dplyr/ggplot2) ----
  "Outcome", "Exposure", "Instrument", "ALLELE0", "ALLELE1", "A1FREQ",
  "beta_exposure", "beta_outcome", "se_exposure", "se_outcome",

  # ---- Derived columns used in plotting / transformations (NSE) ----
  "Estimate", "Lower", "Upper", "Method", "se_ratio",

  # ---- Harmonization / allele manipulation (NSE) ----
  "NEA_exposure", "EA_exposure", "NEA_outcome", "EA_outcome",
  "A1FREQ_exposure", "A1FREQ_outcome",
  "swap_needed", "flip_beta",
  "temp_NEA_outcome", "temp_EA_outcome", "temp_beta_outcome", "temp_A1FREQ_outcome",
  "temp_NEA_exposure", "temp_EA_exposure", "temp_A1FREQ_exposure",
  "temp_NEA_outcome_flip", "temp_EA_outcome_flip", "temp_A1FREQ_outcome_flip",
  "temp_beta_exposure",

  # ---- MR-PRESSO / IVW helpers ----
  "Weights", "mod_noOutliers",

  # ---- JAGS / MR-Horse model symbols (used inside model functions) ----
  "N", "K", "Tx", "R",
  "theta", "mu",
  "by", "sy", "bx", "sx",
  "bx0", "mx0", "vx0", "mx", "sx0",
  "alpha", "phi", "rho", "r", "a", "b", "c", "d", "tau",
  "A", "B", "kappa",
  "inprod", "T"
))
