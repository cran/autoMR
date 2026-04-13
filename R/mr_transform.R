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

#' @noRd
fill_dummy_col <- function(df, col, dummy, label) {
  if (!col %in% names(df)) {
    warning(sprintf(
      "%s column is absent; filling all rows with dummy value %s.",
      label, deparse(dummy)
    ), call. = FALSE)
    df[[col]] <- dummy
  } else {
    na_rows <- is.na(df[[col]])
    n_na    <- sum(na_rows)
    if (n_na > 0) {
      warning(sprintf(
        "%s has %d row(s) with NA; filling with dummy value %s.",
        label, n_na, deparse(dummy)
      ), call. = FALSE)
      df[[col]][na_rows] <- dummy
    }
  }
  df
}

#' @noRd
ensure_dummy_vars <- function(df) {

  # Instrument: auto-generate if column is entirely absent (silent)
  if (!"Instrument" %in% names(df))
    df$Instrument <- paste0("Instrument_", seq_len(nrow(df)))

  # Allele / frequency columns: fill missing or NA rows with dummies, warn once per column
  df <- fill_dummy_col(df, "ALLELE1", "G",  "ALLELE1 (non-effect allele)")
  df <- fill_dummy_col(df, "ALLELE0", "A",  "ALLELE0 (effect allele)")
  df <- fill_dummy_col(df, "A1FREQ",  0.5,  "A1FREQ (effect allele frequency)")

  # Mandatory columns: hard stop if missing
  if (!"Outcome" %in% names(df))
    stop("The data frame is missing the mandatory 'Outcome' column.", call. = FALSE)
  if (!"Exposure" %in% names(df))
    stop("The data frame is missing the mandatory 'Exposure' column.", call. = FALSE)

  return(df)
}

#' @noRd
align_outcome_alleles <- function(df) {
  if (!all(c("NEA_exposure", "EA_exposure", "NEA_outcome", "EA_outcome") %in% names(df)))
    return(df)

  swap_needed <- (df$NEA_exposure != df$NEA_outcome) | (df$EA_exposure != df$EA_outcome)

  old_NEA_outcome <- df$NEA_outcome
  old_EA_outcome  <- df$EA_outcome

  df$NEA_outcome  <- ifelse(swap_needed, old_EA_outcome,  old_NEA_outcome)
  df$EA_outcome   <- ifelse(swap_needed, old_NEA_outcome, old_EA_outcome)
  df$beta_outcome <- ifelse(swap_needed, -df$beta_outcome, df$beta_outcome)

  if ("A1FREQ_outcome" %in% names(df))
    df$A1FREQ_outcome <- ifelse(swap_needed, 1 - df$A1FREQ_outcome, df$A1FREQ_outcome)

  return(df)
}

#' @noRd
apply_beta_sign <- function(df, beta_sign) {
  if (!all(c("NEA_exposure", "EA_exposure") %in% names(df))) return(df)

  # Determine which rows need flipping based on beta_outcome sign
  flip_beta <- if (beta_sign == "positive") {
    df$beta_outcome < 0   # flip rows where outcome beta is negative -> positive
  } else {
    df$beta_outcome > 0   # flip rows where outcome beta is positive -> negative
  }

  # Flip exposure alleles (NEA <-> EA) for flagged rows
  old_NEA_exposure <- df$NEA_exposure
  old_EA_exposure  <- df$EA_exposure
  df$NEA_exposure  <- ifelse(flip_beta, old_EA_exposure,  old_NEA_exposure)
  df$EA_exposure   <- ifelse(flip_beta, old_NEA_exposure, old_EA_exposure)

  # Keep ALLELE1/ALLELE0 in sync immediately
  df$ALLELE1 <- df$NEA_exposure
  df$ALLELE0 <- df$EA_exposure

  # Mirror beta_exposure sign (side effect of strand flip)
  df$beta_exposure <- ifelse(flip_beta, -df$beta_exposure, df$beta_exposure)

  # Apply target sign to beta_outcome: positive -> abs(), negative -> -abs()
  df$beta_outcome <- if (beta_sign == "positive") abs(df$beta_outcome) else -abs(df$beta_outcome)

  # Complement allele frequencies for flagged rows
  if ("A1FREQ_exposure" %in% names(df))
    df$A1FREQ_exposure <- ifelse(flip_beta, 1 - df$A1FREQ_exposure, df$A1FREQ_exposure)
  if ("A1FREQ_outcome" %in% names(df))
    df$A1FREQ_outcome <- ifelse(flip_beta, 1 - df$A1FREQ_outcome, df$A1FREQ_outcome)

  # Flip outcome alleles to match the new exposure strand
  if (all(c("NEA_outcome", "EA_outcome") %in% names(df))) {
    old_NEA_outcome <- df$NEA_outcome
    old_EA_outcome  <- df$EA_outcome
    df$NEA_outcome  <- ifelse(flip_beta, old_EA_outcome,  old_NEA_outcome)
    df$EA_outcome   <- ifelse(flip_beta, old_NEA_outcome, old_EA_outcome)
  }

  return(df)
}

#' @noRd
build_check_df <- function(df) {

  # Resolve SNP identifier (same logic as build_input_df)
  if (!"Instrument" %in% names(df)) {
    if ("SNP" %in% names(df)) {
      df$Instrument <- as.character(df$SNP)
    } else {
      df$Instrument <- paste0("Instrument_", seq_len(nrow(df)))
    }
  }

  # Derive ALLELE1/ALLELE0/A1FREQ from the post-harmonization allele columns
  # so they reflect any flipping that occurred. These are placed at the end
  # of check_df as a preview of what will be passed to input_df.
  if (all(c("NEA_exposure", "EA_exposure") %in% names(df))) {
    df$ALLELE1 <- df$NEA_exposure
    df$ALLELE0 <- df$EA_exposure
  } else {
    if (!"ALLELE1" %in% names(df)) df$ALLELE1 <- "G"
    if (!"ALLELE0" %in% names(df)) df$ALLELE0 <- "A"
  }

  if ("A1FREQ_outcome" %in% names(df)) {
    df$A1FREQ <- df$A1FREQ_outcome
  } else if ("A1FREQ_exposure" %in% names(df)) {
    df$A1FREQ <- df$A1FREQ_exposure
  } else if (!"A1FREQ" %in% names(df)) {
    df$A1FREQ <- 0.5
  }

  # Build the ordered column list:
  #   core columns -> exposure alleles -> outcome alleles -> input_df aliases
  core_cols    <- c("Outcome", "Exposure", "Instrument",
                    "beta_exposure", "se_exposure",
                    "beta_outcome",  "se_outcome")
  exp_cols     <- intersect(c("NEA_exposure", "EA_exposure", "A1FREQ_exposure"), names(df))
  out_cols     <- intersect(c("NEA_outcome",  "EA_outcome",  "A1FREQ_outcome"),  names(df))
  alias_cols   <- c("ALLELE1", "ALLELE0", "A1FREQ")
  # Any extra columns the user had in the original data frame
  extra_cols   <- setdiff(names(df), c(core_cols, exp_cols, out_cols, alias_cols))

  ordered_cols <- c(core_cols, exp_cols, out_cols, extra_cols, alias_cols)
  ordered_cols <- ordered_cols[ordered_cols %in% names(df)]

  return(df[, ordered_cols])
}

#' @noRd
build_input_df <- function(df) {

  # Resolve SNP identifier
  if (!"Instrument" %in% names(df)) {
    if ("SNP" %in% names(df)) {
      df$Instrument <- as.character(df$SNP)
    } else {
      df$Instrument <- paste0("Instrument_", seq_len(nrow(df)))
    }
  }

  # ALLELE1/ALLELE0 from post-harmonization exposure alleles
  if (all(c("NEA_exposure", "EA_exposure") %in% names(df))) {
    df$ALLELE1 <- df$NEA_exposure
    df$ALLELE0 <- df$EA_exposure
  } else {
    if (!"ALLELE1" %in% names(df)) df$ALLELE1 <- "G"
    if (!"ALLELE0" %in% names(df)) df$ALLELE0 <- "A"
  }

  # A1FREQ: prefer outcome, then exposure, then default 0.5
  if ("A1FREQ_outcome" %in% names(df)) {
    df$A1FREQ <- df$A1FREQ_outcome
  } else if ("A1FREQ_exposure" %in% names(df)) {
    df$A1FREQ <- df$A1FREQ_exposure
  } else if (!"A1FREQ" %in% names(df)) {
    df$A1FREQ <- 0.5
  }

  df <- ensure_dummy_vars(df)

  return(df[, c(
    "Outcome", "Exposure", "Instrument",
    "beta_exposure", "se_exposure",
    "beta_outcome",  "se_outcome",
    "ALLELE1", "ALLELE0", "A1FREQ"
  )])
}

# ==============================================================================
# Exported Data Transformation Functions
# ==============================================================================

#' Format vectors into a Mendelian Randomization input data frame
#'
#' Assembles raw vectors into a working data frame, aligns outcome alleles to
#' the exposure strand (Step 1), then standardizes the sign of
#' \code{beta_exposure} (Step 2). Returns both a full working data frame
#' (\code{check_df}) and a slim, renamed input data frame (\code{input_df}),
#' mirroring the output of \code{harmonize_mr_data()}.
#'
#' @param Instrument Character vector of instrument/SNP identifiers.
#' @param beta_exposure Numeric vector of exposure effects.
#' @param se_exposure Numeric vector of exposure standard errors.
#' @param beta_outcome Numeric vector of outcome effects.
#' @param se_outcome Numeric vector of outcome standard errors.
#' @param Outcome Character string or vector for outcome names (Mandatory).
#' @param Exposure Character string or vector for exposure names (Mandatory).
#' @param ALLELE1 Optional character vector for non-effect alleles (NEA) on
#'   the exposure strand (ALLELE1 = NEA_exposure).
#' @param ALLELE0 Optional character vector for effect alleles (EA) on the
#'   exposure strand (ALLELE0 = EA_exposure).
#' @param A1FREQ Optional numeric vector for effect allele frequencies
#'   (exposure dataset).
#' @param ALLELE1_outcome Optional character vector for non-effect alleles in
#'   the outcome dataset. If \code{NULL}, assumed identical to \code{ALLELE1}.
#' @param ALLELE0_outcome Optional character vector for effect alleles in the
#'   outcome dataset. If \code{NULL}, assumed identical to \code{ALLELE0}.
#' @param A1FREQ_outcome Optional numeric vector for effect allele frequencies
#'   in the outcome dataset.
#' @param beta_sign Character string controlling the target sign for
#'   \code{beta_outcome}. One of \code{"positive"} (default, forces
#'   beta_outcome >= 0) or \code{"negative"} (forces beta_outcome <= 0).
#'   When a row is flipped, \code{beta_exposure} is negated as a side effect
#'   and all allele columns are swapped consistently. Ignored when no allele
#'   columns are supplied (strand identity is unknown without allele information).
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{check_df}}{Full working data frame retaining all allele
#'       columns (NEA_exposure, EA_exposure, NEA_outcome, EA_outcome,
#'       A1FREQ_exposure, A1FREQ_outcome) after harmonization. Useful for
#'       quality-checking the harmonization results.}
#'     \item{\code{input_df}}{Slim data frame ready for MR analysis, with
#'       columns: Instrument, beta_exposure, se_exposure, beta_outcome,
#'       se_outcome, Outcome, Exposure, ALLELE1, ALLELE0, A1FREQ.}
#'   }
#' @examples
#' data("fi_49item")
#'
#' # Without allele columns
#' result1 <- format_mr_input(
#'   Instrument    = fi_49item$Instrument,
#'   beta_exposure = fi_49item$beta_exposure,
#'   se_exposure   = fi_49item$se_exposure,
#'   beta_outcome  = fi_49item$beta_outcome,
#'   se_outcome    = fi_49item$se_outcome,
#'   Outcome       = fi_49item$Outcome,
#'   Exposure      = fi_49item$Exposure
#' )
#' head(result1$input_df)
#'
#' # With allele columns (enables alignment + sign standardization)
#' result2 <- format_mr_input(
#'   Instrument      = fi_49item$Instrument,
#'   beta_exposure   = fi_49item$beta_exposure,
#'   se_exposure     = fi_49item$se_exposure,
#'   beta_outcome    = fi_49item$beta_outcome,
#'   se_outcome      = fi_49item$se_outcome,
#'   Outcome         = fi_49item$Outcome,
#'   Exposure        = fi_49item$Exposure,
#'   ALLELE1         = fi_49item$NEA_exposure,
#'   ALLELE0         = fi_49item$EA_exposure,
#'   A1FREQ          = fi_49item$A1FREQ_exposure,
#'   ALLELE1_outcome = fi_49item$NEA_outcome,
#'   ALLELE0_outcome = fi_49item$EA_outcome,
#'   A1FREQ_outcome  = fi_49item$A1FREQ_outcome
#' )
#' head(result2$check_df)
#' head(result2$input_df)
#'
#' # Force all exposure betas to be negative
#' result3 <- format_mr_input(
#'   Instrument    = fi_49item$Instrument,
#'   beta_exposure = fi_49item$beta_exposure,
#'   se_exposure   = fi_49item$se_exposure,
#'   beta_outcome  = fi_49item$beta_outcome,
#'   se_outcome    = fi_49item$se_outcome,
#'   Outcome       = fi_49item$Outcome,
#'   Exposure      = fi_49item$Exposure,
#'   ALLELE1       = fi_49item$NEA_exposure,
#'   ALLELE0       = fi_49item$EA_exposure,
#'   beta_sign     = "negative"
#' )
#' head(result3$input_df)
#' @export
format_mr_input <- function(Instrument, beta_exposure, se_exposure, beta_outcome, se_outcome,
                            Outcome, Exposure,
                            ALLELE1 = NULL, ALLELE0 = NULL, A1FREQ = NULL,
                            ALLELE1_outcome = NULL, ALLELE0_outcome = NULL,
                            A1FREQ_outcome = NULL,
                            beta_sign = c("positive", "negative")) {

  beta_sign <- match.arg(beta_sign)

  df <- data.frame(
    Instrument    = Instrument,
    beta_exposure = beta_exposure,
    se_exposure   = se_exposure,
    beta_outcome  = beta_outcome,
    se_outcome    = se_outcome,
    Outcome       = Outcome,
    Exposure      = Exposure,
    stringsAsFactors = FALSE
  )

  # Map ALLELE1/ALLELE0 -> named exposure allele columns.
  # ALLELE1 = NEA (non-effect allele), ALLELE0 = EA (effect allele).
  if (!is.null(ALLELE1)) { df$NEA_exposure <- ALLELE1; df$ALLELE1 <- ALLELE1 }
  if (!is.null(ALLELE0)) { df$EA_exposure  <- ALLELE0; df$ALLELE0 <- ALLELE0 }
  if (!is.null(A1FREQ))  { df$A1FREQ_exposure <- A1FREQ; df$A1FREQ <- A1FREQ }

  # Populate outcome allele columns. Default to exposure alleles when absent
  # so the alignment step is a no-op (correct behaviour).
  if (!is.null(ALLELE1_outcome)) {
    df$NEA_outcome <- ALLELE1_outcome
  } else if (!is.null(ALLELE1)) {
    df$NEA_outcome <- ALLELE1
  }
  if (!is.null(ALLELE0_outcome)) {
    df$EA_outcome <- ALLELE0_outcome
  } else if (!is.null(ALLELE0)) {
    df$EA_outcome <- ALLELE0
  }
  if (!is.null(A1FREQ_outcome)) df$A1FREQ_outcome <- A1FREQ_outcome

  # Step 1: align outcome alleles to exposure alleles
  df <- align_outcome_alleles(df)

  # Step 2: standardize beta_exposure sign (no-op when allele columns absent)
  df <- apply_beta_sign(df, beta_sign)

  return(list(
    check_df = build_check_df(df),
    input_df = build_input_df(df)
  ))
}

#' Harmonize exposure and outcome SNP data
#'
#' Takes a data frame that already contains all required MR columns, aligns
#' outcome alleles to the exposure strand (Step 1), then standardizes the sign
#' of \code{beta_exposure} (Step 2). Returns both a full working data frame
#' (\code{check_df}) and a slim, renamed input data frame (\code{input_df}),
#' mirroring the output of \code{format_mr_input()}.
#'
#' @param df Data frame containing required columns: \code{beta_exposure},
#'   \code{se_exposure}, \code{beta_outcome}, \code{se_outcome},
#'   \code{Outcome}, \code{Exposure}. Allele columns may be supplied as
#'   \code{NEA_exposure}/\code{EA_exposure}/\code{NEA_outcome}/\code{EA_outcome},
#'   or as \code{ALLELE1} (NEA) / \code{ALLELE0} (EA) which are treated as
#'   exposure alleles.
#' @param Outcome Optional; character string to set as the Outcome name when
#'   the column is absent from \code{df}.
#' @param Exposure Optional; character string to set as the Exposure name when
#'   the column is absent from \code{df}.
#' @param beta_sign Character string controlling the target sign for
#'   \code{beta_outcome} after harmonization. One of \code{"positive"}
#'   (default, forces beta_outcome >= 0) or \code{"negative"} (forces
#'   beta_outcome <= 0). When a row is flipped, \code{beta_exposure} is
#'   negated as a side effect and all allele columns are swapped consistently.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{check_df}}{Full working data frame retaining all allele
#'       columns after harmonization. Useful for quality-checking results.}
#'     \item{\code{input_df}}{Slim data frame ready for MR analysis, with
#'       columns: Instrument, beta_exposure, se_exposure, beta_outcome,
#'       se_outcome, Outcome, Exposure, ALLELE1, ALLELE0, A1FREQ.}
#'   }
#' @import dplyr
#' @examples
#' result1 <- harmonize_mr_data(df = fi_49item)
#' head(result1$check_df)
#' head(result1$input_df)
#'
#' # Force all exposure betas to be negative
#' result2 <- harmonize_mr_data(df = fi_49item, beta_sign = "negative")
#' head(result2$input_df)
#' @export
harmonize_mr_data <- function(df, Outcome = NULL, Exposure = NULL,
                              beta_sign = c("positive", "negative")) {

  beta_sign <- match.arg(beta_sign)

  # Inject metadata if provided via arguments
  if (!is.null(Outcome)  && !"Outcome"  %in% names(df)) df$Outcome  <- Outcome
  if (!is.null(Exposure) && !"Exposure" %in% names(df)) df$Exposure <- Exposure

  # Validate required columns
  req  <- c("beta_exposure", "se_exposure", "beta_outcome", "se_outcome", "Outcome", "Exposure")
  miss <- setdiff(req, names(df))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)

  # Accept ALLELE1/ALLELE0 as aliases for NEA_exposure/EA_exposure
  if (!"NEA_exposure" %in% names(df) && "ALLELE1" %in% names(df))
    df$NEA_exposure <- df$ALLELE1
  if (!"EA_exposure"  %in% names(df) && "ALLELE0" %in% names(df))
    df$EA_exposure  <- df$ALLELE0

  # Default outcome allele columns to exposure alleles when absent
  if (!"NEA_outcome" %in% names(df) && "NEA_exposure" %in% names(df))
    df$NEA_outcome <- df$NEA_exposure
  if (!"EA_outcome"  %in% names(df) && "EA_exposure"  %in% names(df))
    df$EA_outcome  <- df$EA_exposure

  # Step 1: align outcome alleles to exposure alleles
  df <- align_outcome_alleles(df)

  # Step 2: standardize beta_exposure sign (flip strand + mirror outcome)
  df <- apply_beta_sign(df, beta_sign)

  return(list(
    check_df = build_check_df(df),
    input_df = build_input_df(df)
  ))
}
