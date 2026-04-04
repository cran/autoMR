# ==============================================================================
# Forest Plots: Instrument-level and Method-level visualizations
# ==============================================================================

# Internal global variable bindings are now consolidated in R/globals.R

# ==============================================================================
# Internal Helper Functions
# ==============================================================================

#' Extract Method Column Value
#'
#' Helper to fetch a value with fallbacks from the summary dataframe.
#'
#' @param row A row from the summary dataframe.
#' @param method The MR method name.
#' @param what The metric to extract ("Estimate", "Lower", "Upper", or "P").
#' @param eff The effect scale ("Beta", "OR", "HR").
#' @return A numeric value.
#' @noRd
get_method_col <- function(row, method, what, eff) {
  # Columns use clean names (no _Beta suffix); OR/HR values are exponentiated
  # in the data; Lower/Upper/P names are unchanged.
  nm <- switch(what,
               Estimate = if (method == "PRESSO") "Presso" else if (method == "Horse") "Horse" else if (method == "GRIP") "Grip" else method,
               Lower    = if (method == "PRESSO") "Presso_lower" else if (method == "Horse") "Horse_Lower" else if (method == "GRIP") "Grip_Lower" else paste0(method, "_Lower"),
               Upper    = if (method == "PRESSO") "Presso_upper" else if (method == "Horse") "Horse_Upper" else if (method == "GRIP") "Grip_Upper" else paste0(method, "_Upper"),
               P        = if (method == "Egger")  "Egger_P_value" else if (method == "PRESSO") "Presso_p" else if (method == "Horse") "Horse_P" else if (method == "GRIP") "Grip_P" else paste0(method, "_P")
  )

  if (method == "GRIP" && !nm %in% names(row)) {
    alt <- sub("^Grip_", "GRIP_", nm)
    if (alt %in% names(row)) nm <- alt
  }

  if (method == "GRIP" && what %in% c("Lower", "Upper") && !nm %in% names(row)) {
    nm2 <- sub("Grip_Lower$", "Grip_Lower_adj", nm)
    nm2 <- sub("Grip_Upper$", "Grip_Upper_adj", nm2)
    if (nm2 %in% names(row)) nm <- nm2
  }

  if (!nm %in% names(row)) return(NA_real_)
  suppressWarnings(as.numeric(row[[nm]]))
}

# ==============================================================================
# Exported User-Facing Functions
# ==============================================================================

#' Generate Instrument-level Forest Plots for Mendelian Randomization
#'
#' Creates a forest plot for each instrument (SNP) within the MR analysis,
#' including a pooled IVW estimate at the top for comparison.
#'
#' @param MR_input_data Harmonised MR input data frame. Must contain Outcome, Exposure,
#'   Instrument, beta_exposure, se_exposure, beta_outcome, and se_outcome columns.
#' @param report_form Character string or vector indicating the standard
#'   output scale for each outcome (e.g., "Beta", "OR", "HR"). Defaults to "Beta".
#' @param xlim_custom Optional numeric vector of length 2 providing custom
#'   limits for the x-axis. If NULL, limits are determined by the data.
#' @param dot_size Numeric value specifying the size of the points. Default is 2.
#' @param axis_text_size Numeric value specifying the font size for axis labels.
#'   Default is 10.
#' @param axis_title_size Numeric value specifying the font size for axis titles.
#'   Default is 12.
#' @param digits Integer specifying the number of decimal places for labels.
#'   Default is 2.
#' @param plot_width Numeric value specifying the width of exported plot file.
#'   Default is 8.
#' @param plot_height Numeric value specifying the height of exported plot file.
#'   Default is 6.
#' @param label_text_size Numeric value specifying the size of estimate labels
#'   (Beta/OR/HR and 95 percent CI) shown on the plot. Default is 3.
#'
#' @return A \code{GWASForestPlots} object containing one \code{ggplot}
#'   per outcome-exposure pair, with instrument-level (SNP) causal estimates
#'   and a pooled IVW estimate. Use \code{export_forest_plots()} to
#'   write plots to disk with optional filtering.
#' @examples
#' \donttest{
#' data("merged_data")
#' input3 <- harmonize_mr_data(df = merged_data)
#'
#' gwas_plots <- GWAS_forest(
#'   MR_input_data   = input3,
#'   report_form     = c("Beta","OR"),
#'   xlim_custom     = NULL,
#'   dot_size        = 2,
#'   axis_text_size  = 10,
#'   axis_title_size = 12,
#'   digits          = 2,
#'   plot_width      = 8,
#'   plot_height     = 6,
#'   label_text_size = 3
#' )
#'
#' # Retrieve the exact outcome/exposure labels stored in the object
#' out_name <- gwas_plots@outcomes[1]   # "fi_49item"
#' exp_name <- gwas_plots@exposures[1]  # "Zn"
#'
#' # Export all instrument-level plots as PNG (commented — writes to disk)
#' # export_forest_plots(gwas_plots, save_dir = tempdir(), file_type = "jpeg")
#'
#' # Export plots for one outcome only
#' # export_forest_plots(gwas_plots, save_dir = tempdir(), outcome = out_name)
#'
#' # Export plots for one exposure only
#' # export_forest_plots(gwas_plots, save_dir = tempdir(), exposure = exp_name)
#'
#' # Export one specific outcome-exposure pair
#' # export_forest_plots(gwas_plots, save_dir = tempdir(), outcome = out_name, exposure = exp_name)
#' }
#' @import ggplot2
#' @import dplyr
#' @importFrom MendelianRandomization mr_input mr_ivw
#' @export
GWAS_forest <- function(MR_input_data, report_form,
                        xlim_custom = NULL,
                        dot_size = 2,
                        axis_text_size = 10,
                        axis_title_size = 12,
                        digits = 2,
                        plot_width = 8,
                        plot_height = 6,
                        label_text_size = 3) {


  df <- MR_input_data

  unique_outcomes <- unique(df$Outcome)
  if (length(report_form) == 1) report_form <- rep(report_form, length(unique_outcomes))

  outcome_exposure <- df %>% dplyr::distinct(Outcome, Exposure)
  outcome_exposure$effect <- sapply(outcome_exposure$Outcome, function(x) report_form[which(unique_outcomes == x)])

  plots_list    <- list()
  outcomes_vec  <- character()
  exposures_vec <- character()
  for (i in seq_len(nrow(outcome_exposure))) {
    current_outcome  <- outcome_exposure$Outcome[i]
    current_exposure <- outcome_exposure$Exposure[i]
    current_form     <- as.character(outcome_exposure$effect[i])

    df_sub <- df %>%
      dplyr::filter(Outcome == current_outcome, Exposure == current_exposure) %>%
      dplyr::mutate(
        Estimate = beta_outcome / beta_exposure,
        se_ratio = sqrt((se_outcome / beta_exposure)^2 + ((beta_outcome * se_exposure) / (beta_exposure^2))^2),
        Lower    = Estimate - 1.96 * se_ratio,
        Upper    = Estimate + 1.96 * se_ratio
      )

    if (current_form != "Beta") {
      df_sub <- df_sub %>% dplyr::mutate(Estimate = exp(Estimate), Lower = exp(Lower), Upper = exp(Upper))
      ref_line <- 1
    } else {
      ref_line <- 0
    }

    MRInput      <- MendelianRandomization::mr_input(bx = df_sub$beta_exposure, bxse = df_sub$se_exposure, by = df_sub$beta_outcome, byse = df_sub$se_outcome)
    ivw_obj      <- MendelianRandomization::mr_ivw(MRInput)
    pooled_beta  <- if (current_form != "Beta") exp(ivw_obj@Estimate) else ivw_obj@Estimate
    pooled_lower <- if (current_form != "Beta") exp(ivw_obj@CILower) else ivw_obj@CILower
    pooled_upper <- if (current_form != "Beta") exp(ivw_obj@CIUpper) else ivw_obj@CIUpper

    pooled_row <- dplyr::tibble(Instrument = "IVW Pooled", Estimate = pooled_beta, Lower = pooled_lower, Upper = pooled_upper)
    df_plot <- df_sub %>% dplyr::select(Instrument, Estimate, Lower, Upper) %>% dplyr::bind_rows(pooled_row)

    inst_order <- df_plot %>% dplyr::filter(Instrument != "IVW Pooled") %>% dplyr::arrange(Estimate) %>% dplyr::pull(Instrument)
    df_plot <- df_plot %>% dplyr::mutate(Instrument = factor(Instrument, levels = c("IVW Pooled", inst_order)))

    all_vals <- c(df_plot$Estimate, df_plot$Lower, df_plot$Upper)
    data_limits <- if(!is.null(xlim_custom)) xlim_custom else c(min(all_vals, na.rm=TRUE), max(all_vals, na.rm=TRUE))

    fmt_str <- paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]")

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Estimate, y = Instrument)) +
      ggplot2::geom_point(size = dot_size) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lower, xmax = Upper), height = dot_size / 20) +
      ggplot2::geom_vline(xintercept = ref_line, linetype = "dashed") +
      ggplot2::scale_x_continuous(limits = data_limits, expand = ggplot2::expansion(mult = c(0, 0.2))) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf(fmt_str, Estimate, Lower, Upper)), x = data_limits[2], hjust = 0, size = label_text_size) +
      ggplot2::labs(title = paste("Exposure:", current_exposure, "\nOutcome:", current_outcome), x = current_form, y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                     axis.text = ggplot2::element_text(size = axis_text_size), axis.title = ggplot2::element_text(size = axis_title_size),
                     plot.title = ggplot2::element_text(face = "bold", size = axis_title_size + 2))

    plots_list[[paste0(current_outcome, '::', current_exposure)]] <- p
    outcomes_vec  <- c(outcomes_vec,  current_outcome)
    exposures_vec <- c(exposures_vec, current_exposure)
  }
  new_gwas_forest_plots(
    plots     = plots_list,
    outcomes  = outcomes_vec,
    exposures = exposures_vec
  )
}

#' Generate Forest Plots across Multiple MR Methods
#'
#' Creates a forest plot comparing causal estimates across different MR
#' methods (e.g., IVW, RAPS, Egger).
#'
#' @param summary_df MR results data frame, typically the output from run_mr_analysis().
#' @param effect Character string or vector indicating the effect scale
#'   ("Beta", "OR", or "HR").
#' @param custom_xlim Optional numeric vector of length 2 for x-axis limits.
#' @param dot_size Numeric value specifying the point size. Default is 3.
#' @param axis_text_size Numeric value for axis font size.
#' @param axis_title_size Numeric value for title font size.
#' @param pval_text_size Numeric value for p-value label size.
#' @param plot_width Width of exported plot.
#' @param plot_height Height of exported plot.
#' @param clamp_nonpositive Logical; whether non-positive estimates should
#'   be clamped to a small positive value before log-transformation.
#'
#' @return An \code{MRForestPlots} object containing one \code{ggplot}
#'   per outcome-exposure pair, with causal estimates compared across MR
#'   methods. Use \code{export_forest_plots()} to write plots to disk
#'   with optional filtering.
#' @examples
#' \donttest{
#' data("merged_data")
#' input3   <- harmonize_mr_data(df = merged_data)
#' outcome3 <- run_mr_analysis(
#'   MR_input_data     = input3,
#'   outcome.form      = c("Beta","OR"),
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
#'
#' mr_plots <- MR_forest(
#'   summary_df        = outcome3,
#'   effect            = c("Beta","OR"),
#'   custom_xlim       = NULL,
#'   dot_size          = 3,
#'   axis_text_size    = 10,
#'   axis_title_size   = 12,
#'   pval_text_size    = 3,
#'   plot_width        = 8,
#'   plot_height       = 6,
#'   clamp_nonpositive = FALSE
#' )
#'
#' # Retrieve the exact outcome/exposure labels stored in the object
#' out_name <- mr_plots@outcomes[1]   # "fi_49item"
#' exp_name <- mr_plots@exposures[1]  # "Zn"
#'
#' # Export all method-level plots as PDF (commented — writes to disk)
#' export_forest_plots(mr_plots, save_dir = tempdir(), file_type = "png")
#'
#' # Export plots for one outcome only
#' export_forest_plots(mr_plots, save_dir = tempdir(), outcome = out_name)
#'
#' # Export plots for one exposure only
#' export_forest_plots(mr_plots, save_dir = tempdir(), exposure = exp_name)
#'
#' # Export one specific outcome-exposure pair
#' export_forest_plots(mr_plots, save_dir = tempdir(), outcome = out_name, exposure = exp_name)
#' }
#' @import ggplot2
#' @import dplyr
#' @importFrom tidyr drop_na
#' @export
MR_forest <- function(summary_df, effect,
                      custom_xlim = NULL,
                      dot_size = 3,
                      axis_text_size = 10,
                      axis_title_size = 12,
                      pval_text_size = 3,
                      plot_width = 8,
                      plot_height = 6,
                      clamp_nonpositive = FALSE) {


  df <- summary_df

  outcome_exposure <- df %>% dplyr::distinct(Outcome, Exposure)
  unique_outcomes  <- unique(df$Outcome)

  if(length(effect) == 1) effect <- rep(effect, length(unique_outcomes))
  outcome_exposure$effect <- sapply(outcome_exposure$Outcome, function(x) effect[which(unique_outcomes == x)])

  methods       <- c("IVW", "RAPS", "Med", "Egger", "PRESSO", "Horse", "GRIP")
  plots_list    <- list()
  outcomes_vec  <- character()
  exposures_vec <- character()

  for(i in seq_len(nrow(outcome_exposure))) {
    curr_out <- outcome_exposure$Outcome[i]
    curr_exp <- outcome_exposure$Exposure[i]
    curr_eff <- as.character(outcome_exposure$effect[i])
    xlab_val <- if(curr_eff == "Beta") curr_eff else paste0("log(",curr_eff,")")

    df_row <- df %>% dplyr::filter(Outcome == curr_out, Exposure == curr_exp)
    inst_count_label <- df_row$SNPs
    method_df <- data.frame(Method = methods, stringsAsFactors = FALSE)

    est_raw <- sapply(methods, function(m) get_method_col(df_row, m, "Estimate", curr_eff))
    lo_raw  <- sapply(methods, function(m) get_method_col(df_row, m, "Lower",    curr_eff))
    hi_raw  <- sapply(methods, function(m) get_method_col(df_row, m, "Upper",    curr_eff))

    if(curr_eff == "Beta") {
      method_df$Estimate <- est_raw; method_df$Lower <- lo_raw; method_df$Upper <- hi_raw
    } else {
      if (clamp_nonpositive){
        eps <- .Machine$double.eps
        est_raw[!is.na(est_raw) & est_raw <= 0] <- eps
        lo_raw [!is.na(lo_raw) & lo_raw <= 0] <- eps
        hi_raw [!is.na(hi_raw) & hi_raw <= 0] <- eps
      }
      method_df$Estimate <- log(est_raw); method_df$Lower <- log(lo_raw); method_df$Upper <- log(hi_raw)
    }

    method_df$p <- as.numeric(sapply(methods, function(m) get_method_col(df_row, m, "P", curr_eff)))
    method_df$Method <- factor(method_df$Method, levels = methods)

    finite_vals <- unlist(method_df[,2:4]); finite_vals <- finite_vals[is.finite(finite_vals)]
    auto_xlim <- if(length(finite_vals)==0) c(-1,1) else c(min(finite_vals), max(finite_vals))
    xlim_to_use <- if(!is.null(custom_xlim)) custom_xlim else auto_xlim
    range_x <- diff(xlim_to_use); if(range_x <= 0) range_x <- 2
    p_label_x <- xlim_to_use[2] + 0.15 * range_x

    p <- ggplot2::ggplot(method_df, ggplot2::aes(x = Estimate, y = Method)) +
      ggplot2::geom_point(size = dot_size, na.rm = TRUE) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lower, xmax = Upper), height = 0.2, linewidth = dot_size * 0.1, na.rm = TRUE) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::scale_x_continuous(limits = c(xlim_to_use[1], p_label_x + 0.1*range_x)) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(title = paste("Exposure:", curr_exp, "\nOutcome:", curr_out),
                    subtitle = paste("Instruments:", inst_count_label), x = xlab_val, y = "") +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                     axis.text = ggplot2::element_text(size = axis_text_size), axis.title = ggplot2::element_text(size = axis_title_size),
                     plot.title = ggplot2::element_text(face = "bold", size = axis_title_size + 2)) +
      ggplot2::annotate("text", x = p_label_x, y = methods, label = ifelse(is.na(method_df$p), "", paste0("p=", formatC(method_df$p, format="f", digits=3))), hjust = 0, size = pval_text_size)

    plots_list[[paste0(curr_out, '::', curr_exp)]] <- p
    outcomes_vec  <- c(outcomes_vec,  curr_out)
    exposures_vec <- c(exposures_vec, curr_exp)
  }
  new_mr_forest_plots(
    plots     = plots_list,
    outcomes  = outcomes_vec,
    exposures = exposures_vec
  )
}
