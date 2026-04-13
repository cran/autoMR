# ==============================================================================
# Global variables declaration for R CMD check
# ==============================================================================

#' @importFrom utils globalVariables
#' @importFrom methods new setClass setValidity setGeneric setMethod
#' @noRd
.autoMR_imports <- function() NULL

utils::globalVariables(c(
  "Outcome", "Exposure", "Instrument", "ALLELE0", "ALLELE1", "A1FREQ",
  "beta_exposure", "beta_outcome", "se_exposure", "se_outcome"
))

# ==============================================================================
# Internal Helper Functions (Calculation & Lookup Engines)
# ==============================================================================

#' Extract method slope from summary table
#' @noRd
get_slope_from_summary <- function(df_row, m, effect_scale) {
  if (is.null(df_row) || nrow(df_row) == 0) return(NA_real_)

  # Columns are always stored with _Beta suffix regardless of scale;
  # OR/HR values are already exponentiated in the data, so we log them back
  # to get the log-scale slope for plotting.
  col <- if (m == "PRESSO") "Presso"
  else if (m == "Horse") "Horse"
  else if (m == "GRIP") "Grip"
  else m

  if (!col %in% names(df_row)) return(NA_real_)
  val <- suppressWarnings(as.numeric(df_row[[col]]))

  if (effect_scale != "Beta") {
    if (is.na(val) || !is.finite(val) || val <= 0) return(NA_real_)
    return(log(val))
  }
  val
}

#' Extract p-value from summary table
#' @noRd
get_p_from_summary <- function(df_row, m) {
  if (is.null(df_row) || nrow(df_row) == 0) return(NA_real_)
  col <- if (m == "Egger") "Egger_P_value"
  else if (m == "PRESSO") "Presso_p"
  else if (m == "Horse") "Horse_P"
  else if (m == "GRIP") "Grip_P"
  else paste0(m, "_P")
  if (!col %in% names(df_row)) return(NA_real_)
  suppressWarnings(as.numeric(df_row[[col]]))
}

#' Internal single-outcome analysis engine
#' @noRd
valid.output <- function(MR_input_data,
                         outcome.form = "Beta",
                         use_ivw = TRUE,
                         use_raps = TRUE,
                         use_median = TRUE,
                         use_egger = TRUE,
                         use_mr_presso = TRUE,
                         use_mr_horse = TRUE,
                         use_mr_grip  = TRUE,
                         NbDistribution = 1000,
                         SignifThreshold = 0.05,
                         mr_horse_n_iter = 5000,
                         mr_horse_n_burnin = 1000,
                         mr_grip_parameters = NULL) {

  MR_input_data <- ensure_dummy_vars(MR_input_data)
  EXP <- unique(MR_input_data$Exposure)
  results_list <- list()

  for(i in seq_along(EXP)){
    clean_Exposure.i <- as.data.frame(MR_input_data[MR_input_data$Exposure == EXP[i], ])
    exposure_beta <- clean_Exposure.i$beta_exposure
    outcome_beta  <- clean_Exposure.i$beta_outcome
    exposure_se   <- clean_Exposure.i$se_exposure
    outcome_se    <- clean_Exposure.i$se_outcome

    MRInputObject <- MendelianRandomization::mr_input(
      bx = exposure_beta, bxse = exposure_se,
      by = outcome_beta,  byse = outcome_se,
      snps = clean_Exposure.i$Instrument,
      effect_allele = clean_Exposure.i$ALLELE0,
      other_allele = clean_Exposure.i$ALLELE1,
      eaf = clean_Exposure.i$A1FREQ
    )

    result_row <- list()
    result_row[["Outcome"]]  <- unique(MR_input_data$Outcome)[1]
    result_row[["Exposure"]] <- EXP[i]

    # --- IVW ---
    if(use_ivw){
      tryCatch({
        ivw_obj <- MendelianRandomization::mr_ivw(MRInputObject)
        result_row[["Instruments"]] <- ivw_obj@SNPs
        result_row[["IVW"]]  <- ivw_obj@Estimate
        result_row[["IVW_Lower"]] <- ivw_obj@CILower
        result_row[["IVW_Upper"]] <- ivw_obj@CIUpper
        result_row[["IVW_P"]]     <- ivw_obj@Pvalue
        heter <- ivw_obj@Heter.Stat; if(length(heter)==1) heter <- c(heter, NA)
        result_row[["IVW_Q"]]     <- heter[1]
        result_row[["IVW_Q_P"]]   <- heter[2]
        result_row[["Fstat"]]     <- ivw_obj@Fstat
      }, error = function(e){
        for(nm in c("Instruments","IVW","IVW_Lower","IVW_Upper","IVW_P","IVW_Q","IVW_Q_P","Fstat")) result_row[[nm]] <- NA
      })
    } else {
      for(nm in c("Instruments","IVW","IVW_Lower","IVW_Upper","IVW_P","IVW_Q","IVW_Q_P","Fstat")) result_row[[nm]] <- NA
    }

    # --- MR-RAPS ---
    if(use_raps){
      tryCatch({
        r <- mr.raps(exposure_beta, outcome_beta, exposure_se, outcome_se)
        result_row[["RAPS"]]  <- r$beta.hat
        result_row[["RAPS_Lower"]] <- r$beta.hat - stats::qnorm(0.975)*r$beta.se
        result_row[["RAPS_Upper"]] <- r$beta.hat + stats::qnorm(0.975)*r$beta.se
        result_row[["RAPS_P"]]     <- r$beta.p.value
      }, error = function(e){
        for(nm in c("RAPS","RAPS_Lower","RAPS_Upper","RAPS_P")) result_row[[nm]] <- NA
      })
    } else {
      for(nm in c("RAPS","RAPS_Lower","RAPS_Upper","RAPS_P")) result_row[[nm]] <- NA
    }

    if(nrow(clean_Exposure.i) > 2){

      # --- Median ---
      if(use_median){
        tryCatch({
          m <- MendelianRandomization::mr_median(MRInputObject)
          result_row[["Med"]]  <- m@Estimate
          result_row[["Med_Lower"]] <- m@CILower
          result_row[["Med_Upper"]] <- m@CIUpper
          result_row[["Med_P"]]     <- m@Pvalue
        }, error = function(e){
          for(nm in c("Med","Med_Lower","Med_Upper","Med_P")) result_row[[nm]] <- NA
        })
      } else {
        for(nm in c("Med","Med_Lower","Med_Upper","Med_P")) result_row[[nm]] <- NA
      }

      # --- Egger ---
      if(use_egger){
        tryCatch({
          eg <- MendelianRandomization::mr_egger(MRInputObject)
          heter <- eg@Heter.Stat; if(length(heter)==1) heter <- c(heter, NA)
          result_row[["Egger"]]       <- eg@Estimate
          result_row[["Egger_Lower"]]      <- eg@CILower.Est
          result_row[["Egger_Upper"]]      <- eg@CIUpper.Est
          result_row[["Egger_P_value"]]    <- eg@Pvalue.Est
          result_row[["Egger_Q"]]          <- heter[1]
          result_row[["Egger_Q_P"]]        <- heter[2]
          result_row[["I_sq"]]             <- eg@I.sq
          result_row[["Intercept_Est"]]    <- eg@Intercept
          result_row[["Intercept_Lower"]]  <- eg@CILower.Int
          result_row[["Intercept_Upper"]]  <- eg@CIUpper.Int
          result_row[["Intercept_P"]]      <- eg@Pvalue.Int
        }, error = function(e){
          for(nm in c("Egger","Egger_Lower","Egger_Upper","Egger_P_value",
                      "Egger_Q","Egger_Q_P","I_sq","Intercept_Est",
                      "Intercept_Lower","Intercept_Upper","Intercept_P")) result_row[[nm]] <- NA
        })
      } else {
        for(nm in c("Egger","Egger_Lower","Egger_Upper","Egger_P_value",
                    "Egger_Q","Egger_Q_P","I_sq","Intercept_Est",
                    "Intercept_Lower","Intercept_Upper","Intercept_P")) result_row[[nm]] <- NA
      }

      # --- PRESSO ---
      if(use_mr_presso && nrow(clean_Exposure.i) >= 4){
        tryCatch({
          pr <- mr_presso(
            BetaOutcome = "beta_outcome", BetaExposure = "beta_exposure",
            SdOutcome   = "se_outcome",   SdExposure   = "se_exposure",
            data = clean_Exposure.i,
            OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
            NbDistribution = NbDistribution, SignifThreshold = SignifThreshold
          )
          b  <- pr$`Main MR results`[1,3]
          se <- pr$`Main MR results`[1,4]
          p  <- pr$`Main MR results`[1,6]
          result_row[["Presso"]]  <- b
          result_row[["Presso_lower"]] <- b - 1.96*se
          result_row[["Presso_upper"]] <- b + 1.96*se
          result_row[["Presso_p"]]     <- p

          orig <- nrow(clean_Exposure.i)
          if(!is.null(pr$`MR-PRESSO results`$`Outlier Test`) && nrow(pr$`MR-PRESSO results`$`Outlier Test`) > 0){
            nout <- nrow(pr$`MR-PRESSO results`$`Outlier Test`)
            result_row[["Presso_Instruments"]]  <- orig - nout
            result_row[["outlier_Instruments"]] <- paste(row.names(pr$`MR-PRESSO results`$`Outlier Test`), collapse = ",")
          } else {
            result_row[["Presso_Instruments"]]  <- orig
            result_row[["outlier_Instruments"]] <- NA
          }
        }, error = function(e){
          for(nm in c("Presso","Presso_lower","Presso_upper","Presso_p","Presso_Instruments","outlier_Instruments")) result_row[[nm]] <- NA
        })
      } else {
        for(nm in c("Presso","Presso_lower","Presso_upper","Presso_p","Presso_Instruments","outlier_Instruments")) result_row[[nm]] <- NA
      }

      # --- Horse ---
      if(use_mr_horse){
        tryCatch({
          D <- data.frame(betaY = outcome_beta, betaX = exposure_beta,
                          betaYse = outcome_se, betaXse = exposure_se)
          ho <- mr_horse(D, no_ini = 3, variable.names = "theta",
                         n.iter = mr_horse_n_iter, n.burnin = mr_horse_n_burnin)
          b <- as.numeric(ho$MR_Estimate$Estimate[1])
          se <- as.numeric(ho$MR_Estimate$SD[1])
          result_row[["Horse"]]  <- b
          result_row[["Horse_Lower"]] <- b - stats::qnorm(0.975)*se
          result_row[["Horse_Upper"]] <- b + stats::qnorm(0.975)*se
          result_row[["Horse_P"]]     <- 2*(1 - stats::pnorm(abs(b/se)))
        }, error = function(e){
          for(nm in c("Horse","Horse_Lower","Horse_Upper","Horse_P")) result_row[[nm]] <- NA
        })
      } else {
        for(nm in c("Horse","Horse_Lower","Horse_Upper","Horse_P")) result_row[[nm]] <- NA
      }

      # --- GRIP ---
      if(use_mr_grip){
        tryCatch({
          params <- if(is.null(mr_grip_parameters)) default_parameters_mr_grip() else mr_grip_parameters
          gr <- mr_grip(
            b_exp = exposure_beta, b_out = outcome_beta,
            se_exp = exposure_se,  se_out = outcome_se,
            parameters = params
          )
          result_row[["Grip"]]  <- as.numeric(gr$b)
          result_row[["Grip_Lower"]] <- as.numeric(gr$b - 1.96*gr$se)
          result_row[["Grip_Upper"]] <- as.numeric(gr$b + 1.96*gr$se)
          result_row[["Grip_P"]]     <- as.numeric(gr$pval)
        }, error = function(e){
          for(nm in c("Grip","Grip_Lower","Grip_Upper","Grip_P")) result_row[[nm]] <- NA
        })
      } else {
        for(nm in c("Grip","Grip_Lower","Grip_Upper","Grip_P")) result_row[[nm]] <- NA
      }

    } else {
      for(nm in c("Med","Med_Lower","Med_Upper","Med_P",
                  "Egger","Egger_Lower","Egger_Upper","Egger_P_value",
                  "Egger_Q","Egger_Q_P","I_sq","Intercept_Est",
                  "Intercept_Lower","Intercept_Upper","Intercept_P",
                  "Presso","Presso_lower","Presso_upper","Presso_p","Presso_Instruments","outlier_Instruments",
                  "Horse","Horse_Lower","Horse_Upper","Horse_P",
                  "Grip","Grip_Lower","Grip_Upper","Grip_P")) result_row[[nm]] <- NA
    }

    results_list[[i]] <- as.data.frame(result_row, stringsAsFactors = FALSE)
  }

  mr_res <- do.call(rbind, results_list)

  if(all(c("Fstat","IVW_Q_P") %in% colnames(mr_res))){
    mr_res$FLe10 <- if(use_ivw) as.integer(mr_res$Fstat < 10) else NA
    mr_res$SigQ  <- if(use_ivw) as.integer(mr_res$IVW_Q_P < 0.01) else NA
  }
  if("Intercept_P" %in% colnames(mr_res)){
    mr_res$I2Le90     <- if(use_egger) as.integer(mr_res$I_sq < 0.9) else NA
    mr_res$Pleiotropy <- if(use_egger) as.integer(mr_res$Intercept_P < 0.01) else NA
  }

  mr_res$Scale <- outcome.form

  new_order <- c("Outcome","Exposure","Instruments","Scale",
                 "Presso_Instruments","outlier_Instruments","FLe10","SigQ","I2Le90","Pleiotropy",
                 setdiff(colnames(mr_res), c("Outcome","Exposure","Instruments","Scale",
                                             "Presso_Instruments","outlier_Instruments","FLe10","SigQ","I2Le90","Pleiotropy")))
  mr_res <- mr_res[, new_order]

  mr_res
}


#' Internal single-exposure scatter plot data assembler
#'
#' Collects all data and pre-computed slope/p-value parameters needed to
#' draw one scatter plot for a single outcome-exposure pair.  Returns a
#' plain list — no graphics device is opened here.  Rendering is done
#' later by \code{.draw_scatter_plot()}.
#'
#' @noRd
MRplots <- function(MR_input_data,
                    d.title        = NULL,
                    d.subtitle     = NULL,
                    plot.xlab,
                    plot.ylab,
                    outcome_label,
                    exposure_label,
                    methods.plot   = c("IVW","RAPS","Egger","PRESSO","Horse","GRIP"),
                    show.legend    = TRUE,
                    summary_df,
                    effect_scale   = "Beta") {

  method_colors <- c(IVW = "chartreuse3", RAPS = "turquoise3", Egger = "cornflowerblue",
                     PRESSO = "red", Horse = "purple", GRIP = "darkorange")
  method_ltys   <- c(IVW = 3, RAPS = 3, Egger = 4, PRESSO = 5, Horse = 6, GRIP = 2)

  dat    <- subset(MR_input_data, Exposure == exposure_label)
  df_row <- subset(summary_df,   Outcome  == outcome_label & Exposure == exposure_label)

  d.x    <- dat$beta_exposure
  d.y    <- dat$beta_outcome
  d.x.se <- dat$se_exposure
  d.y.se <- dat$se_outcome

  x_ext <- max(abs(d.x + d.x.se), na.rm = TRUE)
  y_ext <- max(abs(d.y + d.y.se), abs(d.y - d.y.se), na.rm = TRUE)

  # Pre-compute per-method slopes and p-values
  slopes  <- list()
  pvalues <- list()
  for (m in methods.plot) {
    slopes[[m]]  <- get_slope_from_summary(df_row, m, effect_scale)
    pvalues[[m]] <- get_p_from_summary(df_row, m)
  }

  # Build y-axis label: append scale annotation for OR/HR
  ylab_scale <- if (effect_scale == "Beta") ""
  else paste0(" (log(", effect_scale, "))")
  slope_label <- if (effect_scale == "Beta") "beta" else paste0("log(", effect_scale, ")")

  list(
    d.x          = d.x,
    d.y          = d.y,
    d.x.se       = d.x.se,
    d.y.se       = d.y.se,
    limx         = c(0, x_ext),
    limy         = c(-y_ext, y_ext),
    xlab         = paste(plot.xlab, exposure_label),
    ylab         = paste0(plot.ylab, ylab_scale, " ", outcome_label),
    title        = if (!is.null(d.title)) d.title else
      paste("Exposure:", exposure_label, "\nOutcome:", outcome_label),
    subtitle     = d.subtitle,
    methods      = methods.plot,
    slopes       = slopes,
    pvalues      = pvalues,
    method_colors = method_colors,
    method_ltys   = method_ltys,
    show.legend  = show.legend,
    slope_label  = slope_label
  )
}

#' Internal scatter plot renderer
#'
#' Draws one scatter plot onto the currently active graphics device using
#' the parameter list produced by \code{MRplots()}.  Must be called while
#' a device is already open.
#'
#' @noRd
.draw_scatter_plot <- function(p) {

  graphics::plot(
    p$d.x, p$d.y,
    xlim = p$limx, ylim = p$limy,
    xlab = p$xlab,
    ylab = p$ylab,
    main = p$title,
    pch  = 16, bty = "L"
  )
  if (!is.null(p$subtitle)) graphics::mtext(p$subtitle)
  graphics::abline(h = 0, lty = 2)
  graphics::abline(v = 0, lty = 2)
  graphics::arrows(p$d.x, p$d.y - p$d.y.se, p$d.x, p$d.y + p$d.y.se,
                   length = 0, angle = 90, code = 3, col = "grey")
  graphics::arrows(p$d.x - p$d.x.se, p$d.y, p$d.x + p$d.x.se, p$d.y,
                   length = 0, angle = 90, code = 3, col = "grey")

  legend_txt  <- character()
  legend_cols <- character()
  legend_ltys <- integer()

  for (m in p$methods) {
    b   <- p$slopes[[m]]
    pv  <- p$pvalues[[m]]
    col <- p$method_colors[[m]]
    lty <- p$method_ltys[[m]]

    if (!is.na(b) && is.finite(b)) {
      graphics::abline(a = 0, b = b, lty = lty, lwd = 2, col = col)
      legend_txt <- c(legend_txt, sprintf("%s: %s=%.3f, p=%.3g", m, p$slope_label, b, pv))
    } else {
      legend_txt <- c(legend_txt, sprintf("%s: failed", m))
    }
    legend_cols <- c(legend_cols, col)
    legend_ltys <- c(legend_ltys, lty)
  }

  if (p$show.legend && length(legend_txt) > 0) {
    graphics::legend("topright", legend = legend_txt,
                     col = legend_cols, lty = legend_ltys,
                     lwd = 2, bg = "white", cex = 0.9)
  }

  invisible(NULL)
}
# ==============================================================================
# Exported User-Facing Functions
# ==============================================================================

#' Run MR Analysis for Multiple Outcomes
#'
#' Performs causal inference analysis using multiple Mendelian Randomization
#' (MR) methods across one or more outcomes and exposures. Returns a combined
#' results data frame. To save the output, use standard R functions such as
#' \code{write.csv()} or \code{saveRDS()} on the returned object.
#'
#' @param MR_input_data Harmonised MR input data frame. Must contain Outcome
#'   and Exposure columns.
#' @param outcome.form Character vector indicating the effect scale for each
#'   outcome: \code{"Beta"}, \code{"OR"} (odds ratio), or \code{"HR"}
#'   (hazard ratio). A single value is recycled across all outcomes.
#'   Defaults to \code{"Beta"}.
#' @param use_ivw Logical; whether to run the Inverse Variance Weighted
#'   (IVW) method. Default is \code{TRUE}.
#' @param use_raps Logical; whether to run the Robust Adjusted Profile Score
#'   (MR-RAPS) method. Default is \code{TRUE}.
#' @param use_median Logical; whether to run the Weighted Median method.
#'   Default is \code{TRUE}.
#' @param use_egger Logical; whether to run MR-Egger regression.
#'   Default is \code{TRUE}.
#' @param use_mr_presso Logical; whether to run the Mendelian Randomization
#'   Pleiotropy RESidual Sum and Outlier (MR-PRESSO) method.
#'   Default is \code{TRUE}.
#' @param use_mr_horse Logical; whether to run the MR-Horse method.
#'   Default is \code{TRUE}.
#' @param use_mr_grip Logical; whether to run the Generalized Regression with
#'   Instrument Pairs (MR-GRIP) method. Default is \code{TRUE}.
#' @param NbDistribution Integer; number of simulated distributions for
#'   MR-PRESSO. Default is \code{1000}.
#' @param SignifThreshold Numeric; significance threshold for the MR-PRESSO
#'   outlier test. Default is \code{0.05}.
#' @param mr_horse_n_iter Integer; number of Markov chain Monte Carlo (MCMC)
#'   iterations for MR-Horse. Default is \code{5000}.
#' @param mr_horse_n_burnin Integer; number of MCMC burn-in samples for
#'   MR-Horse. Default is \code{1000}.
#' @param mr_grip_parameters List of additional parameters passed to MR-GRIP.
#'   If \code{NULL}, default parameters are used.
#'
#' @return A data frame combining results across all outcomes and exposures.
#'   Each row represents one outcome-exposure pair. Columns include estimates,
#'   confidence intervals (CI), and p-values for each method, together with
#'   diagnostic flags (e.g., F-statistic below 10, significant heterogeneity).
#'   Use \code{write.csv()} or \code{saveRDS()} to save the returned object.
#'
#' @examples
#' data("fi_49item")
#' input1 <- harmonize_mr_data(df = fi_49item)$input_df
#' outcome1 <- run_mr_analysis(
#'   MR_input_data     = input1,
#'   outcome.form      = "Beta",
#'   use_ivw           = TRUE,
#'   use_raps          = FALSE,
#'   use_median        = FALSE,
#'   use_egger         = FALSE,
#'   use_mr_presso     = FALSE,
#'   use_mr_horse      = FALSE,
#'   use_mr_grip       = FALSE,
#'   NbDistribution    = 1000,
#'   SignifThreshold   = 0.05,
#'   mr_horse_n_iter   = 5000,
#'   mr_horse_n_burnin = 1000,
#'   mr_grip_parameters = NULL
#' )
#'
#' \donttest{
#' data("fried_frailty")
#' input2 <- harmonize_mr_data(df = fried_frailty)$input_df
#' outcome2 <- run_mr_analysis(
#'   MR_input_data     = input2,
#'   outcome.form      = "OR",
#'   use_ivw           = TRUE,
#'   use_raps          = TRUE,
#'   use_median        = TRUE,
#'   use_egger         = TRUE,
#'   use_mr_presso     = TRUE,
#'   use_mr_horse      = TRUE,
#'   use_mr_grip       = TRUE,
#'   NbDistribution    = 1000,
#'   SignifThreshold   = 0.05,
#'   mr_horse_n_iter   = 5000,
#'   mr_horse_n_burnin = 1000,
#'   mr_grip_parameters = NULL
#' )
#' }
#'
#' \donttest{
#' data("merged_data")
#' input3 <- harmonize_mr_data(df = merged_data)$input_df
#' outcome3 <- run_mr_analysis(
#'   MR_input_data     = input3,
#'   outcome.form      = c("Beta","OR"), ## First outcome use Beta and second outcome use OR
#'   use_ivw           = TRUE,
#'   use_raps          = TRUE,
#'   use_median        = TRUE,
#'   use_egger         = TRUE,
#'   use_mr_presso     = TRUE,
#'   use_mr_horse      = TRUE,
#'   use_mr_grip       = TRUE,
#'   NbDistribution    = 1000,
#'   SignifThreshold   = 0.05,
#'   mr_horse_n_iter   = 5000,
#'   mr_horse_n_burnin = 1000,
#'   mr_grip_parameters = NULL
#' )
#' }
#' @export
run_mr_analysis <- function(MR_input_data,
                            outcome.form      = NULL,
                            use_ivw           = TRUE,
                            use_raps          = TRUE,
                            use_median        = TRUE,
                            use_egger         = TRUE,
                            use_mr_presso     = TRUE,
                            use_mr_horse      = TRUE,
                            use_mr_grip       = TRUE,
                            NbDistribution    = 1000,
                            SignifThreshold   = 0.05,
                            mr_horse_n_iter   = 5000,
                            mr_horse_n_burnin = 1000,
                            mr_grip_parameters = NULL) {

  outcomes <- unique(MR_input_data$Outcome)
  if (is.null(outcome.form)) outcome.form <- rep("Beta", length(outcomes))
  if (length(outcome.form) == 1) outcome.form <- rep(outcome.form, length(outcomes))

  results_list <- vector("list", length(outcomes))
  for (i in seq_along(outcomes)) {
    Outcome.i <- as.data.frame(MR_input_data[MR_input_data$Outcome == outcomes[i], ])
    results_list[[i]] <- valid.output(
      MR_input_data     = Outcome.i,
      outcome.form      = outcome.form[i],
      use_ivw           = use_ivw,
      use_raps          = use_raps,
      use_median        = use_median,
      use_egger         = use_egger,
      use_mr_presso     = use_mr_presso,
      use_mr_horse      = use_mr_horse,
      use_mr_grip       = use_mr_grip,
      NbDistribution    = NbDistribution,
      SignifThreshold   = SignifThreshold,
      mr_horse_n_iter   = mr_horse_n_iter,
      mr_horse_n_burnin = mr_horse_n_burnin,
      mr_grip_parameters = mr_grip_parameters
    )
  }

  mr_combined <- do.call(rbind, results_list)

  # Apply OR/HR exponentiation and column renaming per-outcome after rbind
  exp_cols_beta  <- c("IVW","RAPS","Med","Egger","Presso","Horse","Grip",
                      "IVW_Lower","RAPS_Lower","Med_Lower","Egger_Lower","Presso_lower","Horse_Lower","Grip_Lower",
                      "IVW_Upper","RAPS_Upper","Med_Upper","Egger_Upper","Presso_upper","Horse_Upper","Grip_Upper")

  for (i in seq_along(outcomes)) {
    form <- outcome.form[i]
    if (form %in% c("OR", "HR")) {
      rows <- mr_combined$Outcome == outcomes[i]
      for (col in exp_cols_beta) {
        if (col %in% colnames(mr_combined))
          mr_combined[rows, col] <- exp(mr_combined[rows, col])
      }
    }
  }

  # Rename _Beta columns to _OR or _HR where needed.
  # Since different outcomes may use different scales, we keep _Beta as the
  mr_combined
}

#' Plot MR Scatter Plots for Multiple Outcomes and Exposures
#'
#' Generates one scatter plot per outcome-exposure pair using base R
#' graphics, with a regression line overlaid for each requested
#' Mendelian Randomization (MR) method. Plot parameters are stored in an
#' \code{MRScatterPlots} S4 object and rendered on demand at export
#' time, so no files are written to disk during this call.
#'
#' @param MR_input_data Harmonised MR input data frame. Must contain
#'   \code{Outcome} and \code{Exposure} columns.
#' @param plot.xlab Character string; prefix for the x-axis label.
#'   Default is \code{"Exposure"}.
#' @param plot.ylab Character string; prefix for the y-axis label.
#'   Default is \code{"Outcome"}.
#' @param methods.plot Character vector of MR methods to overlay as
#'   regression lines. Supported values: \code{"IVW"}, \code{"RAPS"},
#'   \code{"Egger"}, \code{"PRESSO"}, \code{"Horse"}, \code{"GRIP"}.
#' @param NbDistribution_presso Integer; number of simulated distributions
#'   for on-the-fly MR-PRESSO calculation. Default is \code{1000}.
#' @param SignifThreshold_presso Numeric; significance threshold for
#'   on-the-fly MR-PRESSO outlier test. Default is \code{0.05}.
#' @param mr_horse_n_iter Integer; number of Markov chain Monte Carlo (MCMC)
#'   iterations for on-the-fly MR-Horse. Default is \code{5000}.
#' @param mr_horse_n_burnin Integer; number of MCMC burn-in samples for
#'   on-the-fly MR-Horse. Default is \code{1000}.
#' @param show.legend Logical; whether to annotate each plot with method
#'   labels, beta estimates, and p-values. Default is \code{TRUE}.
#' @param summary_df Optional data frame of pre-calculated results from
#'   \code{run_mr_analysis()}. When supplied together with
#'   \code{use_df_results = TRUE}, avoids re-running the analysis.
#' @param effect_scale Character string matching the scale used in
#'   \code{summary_df}: \code{"Beta"}, \code{"OR"}, or \code{"HR"}.
#'   Default is \code{"Beta"}.
#' @param use_df_results Logical; if \code{TRUE} and \code{summary_df} is
#'   provided, method slopes are read from \code{summary_df} instead of
#'   being re-calculated. Default is \code{TRUE}.
#'
#' @return An \code{MRScatterPlots} object containing one
#'   \code{recordedplot} per outcome-exposure pair, together with outcome
#'   and exposure metadata. Use \code{export_scatter_plots()} to write
#'   plots to disk with optional filtering by outcome, exposure, or both.
#'
#' @examples
#' data("merged_data")
#' input3  <- harmonize_mr_data(df = merged_data)$input_df
#' outcome3 <- run_mr_analysis(
#'   MR_input_data     = input3,
#'   outcome.form      = c("Beta","OR"),
#'   use_ivw           = TRUE,
#'   use_raps          = FALSE,
#'   use_median        = FALSE,
#'   use_egger         = FALSE,
#'   use_mr_presso     = FALSE,
#'   use_mr_horse      = FALSE,
#'   use_mr_grip       = FALSE,
#'   NbDistribution    = 1000,
#'   SignifThreshold   = 0.05,
#'   mr_horse_n_iter   = 5000,
#'   mr_horse_n_burnin = 1000,
#'   mr_grip_parameters = NULL
#' )
#' \donttest{
#' # Pass pre-calculated results to avoid rerunning the analysis
#' plots <- plot_mr_scatter(
#'   MR_input_data  = input3,
#'   summary_df     = outcome3,
#'   use_df_results = TRUE
#' )
#'
#' # Inspect the object; for fi_49item this prints:
#' #   [1] fi_49item :: Zn
#'
#' # Retrieve the exact outcome/exposure labels stored in the object
#' out_name <- plots@outcomes[1]   # "fi_49item"
#' exp_name <- plots@exposures[1]  # "Zn"
#'
#' # Export all plots as PDF (commented — writes to disk)
#' export_scatter_plots(plots, save_dir = tempdir(), file_type = "pdf")
#'
#' # Export one outcome only
#' export_scatter_plots(plots, save_dir = tempdir(), outcome = out_name)
#'
#' # Export one exposure only
#' export_scatter_plots(plots, save_dir = tempdir(), exposure = exp_name)
#'
#' # Export one specific pair
#' # export_scatter_plots(plots, save_dir = tempdir(), outcome = out_name, exposure = exp_name)
#'
#' # Export as PNG instead
#' export_scatter_plots(plots, save_dir = tempdir(), file_type = "png")
#'
#' }
#' @importFrom graphics plot abline arrows legend mtext
#' @export
plot_mr_scatter <- function(MR_input_data,
                            plot.xlab             = "Exposure",
                            plot.ylab             = "Outcome",
                            methods.plot          = c("IVW","RAPS","Egger","PRESSO","Horse","GRIP"),
                            NbDistribution_presso = 1000,
                            SignifThreshold_presso = 0.05,
                            mr_horse_n_iter       = 5000,
                            mr_horse_n_burnin     = 1000,
                            show.legend           = TRUE,
                            summary_df            = NULL,
                            effect_scale          = "Beta",
                            use_df_results        = TRUE) {

  MR_input_data <- ensure_dummy_vars(MR_input_data)

  # Compute summary results once if not supplied or not requested
  if (is.null(summary_df) || !use_df_results) {
    summary_df <- valid.output(
      MR_input_data,
      outcome.form   = effect_scale,
      NbDistribution = NbDistribution_presso,
      SignifThreshold = SignifThreshold_presso,
      mr_horse_n_iter   = mr_horse_n_iter,
      mr_horse_n_burnin = mr_horse_n_burnin
    )
  }

  plots_list <- list()
  outcomes_vec  <- character()
  exposures_vec <- character()

  for (out in unique(MR_input_data$Outcome)) {
    sub_data <- subset(MR_input_data, Outcome == out)
    # Use per-outcome Scale from summary_df if available, else fall back to effect_scale
    out_scale <- if (!is.null(summary_df) && "Scale" %in% colnames(summary_df)) {
      sc <- summary_df$Scale[summary_df$Outcome == out]
      if (length(sc) > 0 && !is.na(sc[1])) sc[1] else effect_scale
    } else effect_scale
    for (ex in unique(sub_data$Exposure)) {
      key <- paste0(out, "::", ex)
      plots_list[[key]] <- MRplots(
        MR_input_data  = sub_data,
        outcome_label  = out,
        exposure_label = ex,
        plot.xlab      = plot.xlab,
        plot.ylab      = plot.ylab,
        methods.plot   = methods.plot,
        show.legend    = show.legend,
        summary_df     = summary_df,
        effect_scale   = out_scale
      )
      outcomes_vec  <- c(outcomes_vec,  out)
      exposures_vec <- c(exposures_vec, ex)
    }
  }

  new_mr_scatter_plots(
    plots     = plots_list,
    outcomes  = outcomes_vec,
    exposures = exposures_vec
  )
}
