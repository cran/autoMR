# ==============================================================================
# Global variables declaration for R CMD check
# ==============================================================================
utils::globalVariables(c(
  "Instrument", "NEA_exposure", "EA_exposure", "NEA_outcome", "EA_outcome",
  "A1FREQ_outcome", "swap_needed", "flip_beta", "temp_NEA_outcome",
  "temp_EA_outcome", "temp_beta_outcome", "temp_A1FREQ_outcome",
  "temp_NEA_exposure", "temp_EA_exposure", "temp_A1FREQ_exposure",
  "temp_NEA_outcome_flip", "temp_EA_outcome_flip", "temp_A1FREQ_outcome_flip",
  "temp_beta_exposure", "A1FREQ_exposure", "beta_exposure", "beta_outcome"
))

# ==============================================================================
# Internal Helper Functions
# ==============================================================================

#' Ensure required columns exist, creating dummies only for optional vars
#' @noRd
ensure_dummy_vars <- function(df) {
  # Optional columns: auto-generate if missing
  if (!"Instrument" %in% names(df)) df$Instrument <- paste0("Instrument_", seq_len(nrow(df)))
  if (!"A1FREQ" %in% names(df)) df$A1FREQ <- 0.5
  if (!"ALLELE0" %in% names(df)) df$ALLELE0 <- "A"
  if (!"ALLELE1" %in% names(df)) df$ALLELE1 <- "G"

  # Mandatory columns: stop if missing
  if (!"Outcome" %in% names(df)) {
    stop("The data frame is missing the mandatory 'Outcome' column.", call. = FALSE)
  }
  if (!"Exposure" %in% names(df)) {
    stop("The data frame is missing the mandatory 'Exposure' column.", call. = FALSE)
  }

  return(df)
}

# ==============================================================================
# Exported Data Transformation Functions
# ==============================================================================

#' Format vectors into a Mendelian Randomization input data frame
#'
#' @param Instrument Character vector of instrument/SNP identifiers.
#' @param beta_exposure Numeric vector of exposure effects.
#' @param se_exposure Numeric vector of exposure standard errors.
#' @param beta_outcome Numeric vector of outcome effects.
#' @param se_outcome Numeric vector of outcome standard errors.
#' @param Outcome Character string or vector for outcome names (Mandatory).
#' @param Exposure Character string or vector for exposure names (Mandatory).
#' @param ALLELE1 Optional character vector for non-effect alleles.
#' @param ALLELE0 Optional character vector for effect alleles.
#' @param A1FREQ Optional numeric vector for effect allele frequencies.
#'
#' @return A data frame formatted for MR analysis.
#' @examples
#' data("fi_49item")
#'
#' input_test1 <- format_mr_input(
#'   Instrument = fi_49item$Instrument,
#'   beta_exposure = fi_49item$beta_exposure,
#'   se_exposure = fi_49item$se_exposure,
#'   beta_outcome = fi_49item$beta_outcome,
#'   se_outcome = fi_49item$se_outcome,
#'   Outcome = fi_49item$Outcome,
#'   Exposure = fi_49item$Exposure
#' )
#' head(input_test1)
#'
#' input_test2 <- format_mr_input(
#'   Instrument = fi_49item$Instrument,
#'   beta_exposure = fi_49item$beta_exposure,
#'   se_exposure = fi_49item$se_exposure,
#'   beta_outcome = fi_49item$beta_outcome,
#'   se_outcome = fi_49item$se_outcome,
#'   Outcome = fi_49item$Outcome,
#'   Exposure = fi_49item$Exposure,
#'   ALLELE1 = fi_49item$NEA_outcome,
#'   ALLELE0 = fi_49item$EA_outcome,
#'   A1FREQ = fi_49item$A1FREQ_outcome
#' )
#' head(input_test2)
#' @export
format_mr_input <- function(Instrument, beta_exposure, se_exposure, beta_outcome, se_outcome,
                            Outcome, Exposure,
                            ALLELE1 = NULL, ALLELE0 = NULL, A1FREQ = NULL) {

  df <- data.frame(
    Instrument = Instrument,
    beta_exposure = beta_exposure,
    se_exposure = se_exposure,
    beta_outcome = beta_outcome,
    se_outcome = se_outcome,
    Outcome = Outcome,
    Exposure = Exposure,
    stringsAsFactors = FALSE
  )

  if(!is.null(ALLELE1)) df$ALLELE1 <- ALLELE1
  if(!is.null(ALLELE0)) df$ALLELE0 <- ALLELE0
  if(!is.null(A1FREQ)) df$A1FREQ <- A1FREQ

  return(ensure_dummy_vars(df))
}

#' Harmonize exposure and outcome SNP data
#'
#' Harmonizes alleles and ensures consistent effect direction between exposure
#' and outcome datasets. Standardizes exposure effects to be positive.
#'
#' @param df Data frame containing required columns (beta_outcome, se_outcome,
#' beta_exposure, se_exposure) and metadata (Outcome, Exposure).
#' @param Outcome Optional; Character string to set as the Outcome name.
#' @param Exposure Optional; Character string to set as the Exposure name.
#' @return A harmonized data frame with 'Instrument' and allele columns.
#' @import dplyr
#' @examples
#' input_test3 <- harmonize_mr_data(df = fi_49item)
#' head(input_test3)
#' @export
harmonize_mr_data <- function(df, Outcome = NULL, Exposure = NULL) {

  # Inject metadata if provided via arguments
  if (!is.null(Outcome)  && !"Outcome"  %in% names(df)) df$Outcome  <- Outcome
  if (!is.null(Exposure) && !"Exposure" %in% names(df)) df$Exposure <- Exposure

  # Required numeric columns
  req <- c("beta_exposure", "se_exposure", "beta_outcome", "se_outcome", "Outcome", "Exposure")
  miss <- setdiff(req, names(df))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)

  # ---- 1) Align outcome alleles to exposure alleles (swap outcome if mismatched) ----
  if (all(c("NEA_exposure", "EA_exposure", "NEA_outcome", "EA_outcome") %in% names(df))) {

    swap_needed <- (df$NEA_exposure != df$NEA_outcome) | (df$EA_exposure != df$EA_outcome)

    old_NEA_outcome <- df$NEA_outcome
    old_EA_outcome  <- df$EA_outcome

    df$NEA_outcome  <- ifelse(swap_needed, old_EA_outcome, old_NEA_outcome)
    df$EA_outcome   <- ifelse(swap_needed, old_NEA_outcome, old_EA_outcome)
    df$beta_outcome <- ifelse(swap_needed, -df$beta_outcome, df$beta_outcome)

    if ("A1FREQ_outcome" %in% names(df)) {
      df$A1FREQ_outcome <- ifelse(swap_needed, 1 - df$A1FREQ_outcome, df$A1FREQ_outcome)
    }
  }

  # ---- 2) Force exposure beta >= 0 (flip exposure; flip outcome sign accordingly) ----
  if (all(c("NEA_exposure", "EA_exposure") %in% names(df))) {
    flip_beta <- df$beta_exposure < 0

    old_NEA_exposure <- df$NEA_exposure
    old_EA_exposure  <- df$EA_exposure

    df$NEA_exposure  <- ifelse(flip_beta, old_EA_exposure, old_NEA_exposure)
    df$EA_exposure   <- ifelse(flip_beta, old_NEA_exposure, old_EA_exposure)

    df$beta_exposure <- abs(df$beta_exposure)
    df$beta_outcome  <- ifelse(flip_beta, -df$beta_outcome, df$beta_outcome)

    if ("A1FREQ_exposure" %in% names(df)) {
      df$A1FREQ_exposure <- ifelse(flip_beta, 1 - df$A1FREQ_exposure, df$A1FREQ_exposure)
    }
    if ("A1FREQ_outcome" %in% names(df)) {
      df$A1FREQ_outcome <- ifelse(flip_beta, 1 - df$A1FREQ_outcome, df$A1FREQ_outcome)
    }

    # If outcome allele columns exist, keep them consistent with exposure after flip
    if (all(c("NEA_outcome", "EA_outcome") %in% names(df))) {
      old_NEA_outcome2 <- df$NEA_outcome
      old_EA_outcome2  <- df$EA_outcome
      df$NEA_outcome   <- ifelse(flip_beta, old_EA_outcome2, old_NEA_outcome2)
      df$EA_outcome    <- ifelse(flip_beta, old_NEA_outcome2, old_EA_outcome2)
    }
  }

  # ---- 3) Build final standardized columns ----

  # Instrument
  if ("SNP" %in% names(df)) {
    df$Instrument <- as.character(df$SNP)
  } else if (!"Instrument" %in% names(df)) {
    df$Instrument <- paste0("Instrument_", seq_len(nrow(df)))
  }

  # Alleles -> ALLELE1/ALLELE0
  if (all(c("NEA_exposure", "EA_exposure") %in% names(df))) {
    df$ALLELE1 <- df$NEA_exposure
    df$ALLELE0 <- df$EA_exposure
  } else {
    if (!"ALLELE1" %in% names(df)) df$ALLELE1 <- "G"
    if (!"ALLELE0" %in% names(df)) df$ALLELE0 <- "A"
  }

  # A1FREQ (prefer outcome freq if present; else A1FREQ; else 0.5)
  if ("A1FREQ_outcome" %in% names(df)) {
    df$A1FREQ <- df$A1FREQ_outcome
  } else if (!"A1FREQ" %in% names(df)) {
    df$A1FREQ <- 0.5
  }

  # Return only requested columns in requested order
  res <- df[, c(
    "Instrument",
    "beta_exposure", "se_exposure",
    "beta_outcome",  "se_outcome",
    "Outcome", "Exposure",
    "ALLELE1", "ALLELE0", "A1FREQ"
  )]

  return(res)
}
