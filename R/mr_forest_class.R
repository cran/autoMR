# ==============================================================================
# S4 Classes: GWASForestPlots and MRForestPlots
#
# Typed containers for the ggplot objects produced by GWAS_forest() and
# MR_forest() respectively. Each class stores plots alongside outcome and
# exposure metadata so that filtering and export work on structured fields
# rather than string-parsing of list names.
#
# Both classes share the same slot structure and export interface, but are
# export() can generate the correct filename prefixes automatically.
# ==============================================================================

# ------------------------------------------------------------------------------
# GWASForestPlots — instrument-level forest plots from GWAS_forest()
# ------------------------------------------------------------------------------

#' S4 class to store instrument-level MR forest plots
#'
#' A \code{GWASForestPlots} object is returned by \code{GWAS_forest()}.
#' It holds one \code{ggplot} object per outcome-exposure pair — each plot
#' shows per-instrument (SNP-level) causal estimates alongside a pooled IVW
#' estimate — together with the metadata needed to filter and export subsets.
#'
#' @slot plots     Named list of \code{ggplot} objects. Names follow the
#'   pattern \code{"Outcome::Exposure"}.
#' @slot outcomes  Character vector of outcome labels, parallel to
#'   \code{plots}.
#' @slot exposures Character vector of exposure labels, parallel to
#'   \code{plots}.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{export_forest_plots(object, save_dir, file_type, width,
#'     height, dpi, outcome, exposure)}}{
#'     Saves plots to disk. Optionally filter by \code{outcome},
#'     \code{exposure}, or both to export a subset.}
#' }
#' @importFrom methods setClass setValidity new
#' @export
setClass("GWASForestPlots", representation(
  plots     = "list",
  outcomes  = "character",
  exposures = "character"
))

setValidity("GWASForestPlots", function(object) {
  n <- length(object@plots)
  if (length(object@outcomes)  != n)
    return("'outcomes' must have the same length as 'plots'.")
  if (length(object@exposures) != n)
    return("'exposures' must have the same length as 'plots'.")
  if (n > 0 && is.null(names(object@plots)))
    return("'plots' must be a named list.")
  TRUE
})

# ------------------------------------------------------------------------------
# MRForestPlots — method-level forest plots from MR_forest()
# ------------------------------------------------------------------------------

#' S4 class to store method-level MR forest plots
#'
#' An \code{MRForestPlots} object is returned by \code{MR_forest()}.
#' It holds one \code{ggplot} object per outcome-exposure pair — each plot
#' compares causal estimates across MR methods (IVW, RAPS, Egger, etc.) —
#' together with the metadata needed to filter and export subsets.
#'
#' @slot plots     Named list of \code{ggplot} objects. Names follow the
#'   pattern \code{"Outcome::Exposure"}.
#' @slot outcomes  Character vector of outcome labels, parallel to
#'   \code{plots}.
#' @slot exposures Character vector of exposure labels, parallel to
#'   \code{plots}.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{export_forest_plots(object, save_dir, file_type, width,
#'     height, dpi, outcome, exposure)}}{
#'     Saves plots to disk. Optionally filter by \code{outcome},
#'     \code{exposure}, or both to export a subset.}
#' }
#' @export
setClass("MRForestPlots", representation(
  plots     = "list",
  outcomes  = "character",
  exposures = "character"
))

setValidity("MRForestPlots", function(object) {
  n <- length(object@plots)
  if (length(object@outcomes)  != n)
    return("'outcomes' must have the same length as 'plots'.")
  if (length(object@exposures) != n)
    return("'exposures' must have the same length as 'plots'.")
  if (n > 0 && is.null(names(object@plots)))
    return("'plots' must be a named list.")
  TRUE
})

# ------------------------------------------------------------------------------
# Internal constructors (not exported)
# ------------------------------------------------------------------------------

#' @noRd
new_gwas_forest_plots <- function(plots, outcomes, exposures) {
  new("GWASForestPlots",
      plots     = plots,
      outcomes  = outcomes,
      exposures = exposures)
}

#' @noRd
new_mr_forest_plots <- function(plots, outcomes, exposures) {
  new("MRForestPlots",
      plots     = plots,
      outcomes  = outcomes,
      exposures = exposures)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# export_forest_plots generic and methods
# ------------------------------------------------------------------------------

#' Export MR Forest Plots to Disk
#'
#' Saves plots stored in a \code{GWASForestPlots} or
#' \code{MRForestPlots} object to a directory using
#' \code{ggplot2::ggsave()}. You can export all plots or filter to a specific
#' outcome, a specific exposure, or a specific outcome-exposure pair.
#'
#' @param object    A \code{GWASForestPlots} or \code{MRForestPlots} object
#'   returned by \code{GWAS_forest()} or \code{MR_forest()}.
#' @param save_dir  Character string; directory to write files into.
#'   Must already exist. Defaults to \code{tempdir()}.
#' @param file_type Character string; output format passed to
#'   \code{ggplot2::ggsave()}. One of \code{"pdf"}, \code{"png"},
#'   \code{"jpeg"}, or \code{"tiff"}. Default is \code{"pdf"}.
#' @param width     Numeric; plot width in inches. Default is \code{8}.
#' @param height    Numeric; plot height in inches. Default is \code{6}.
#' @param dpi       Integer; resolution in dots per inch (ignored for PDF).
#'   Default is \code{300}.
#' @param outcome   Optional character string; if supplied, only plots
#'   whose outcome matches this value are exported.
#' @param exposure  Optional character string; if supplied, only plots
#'   whose exposure matches this value are exported.
#'
#' @details Filtering logic:
#' \itemize{
#'   \item \strong{Both \code{outcome} and \code{exposure} supplied} — exports
#'     the single plot for that exact pair.
#'   \item \strong{\code{outcome} only} — exports all exposures for that
#'     outcome.
#'   \item \strong{\code{exposure} only} — exports all outcomes for that
#'     exposure.
#'   \item \strong{Neither supplied} — exports every plot in the object.
#' }
#' File names are prefixed with \code{Instrument_forest_} for
#' \code{GWASForestPlots} objects and \code{Method_forest_} for
#' \code{MRForestPlots} objects, followed by
#' \code{<Outcome>_<Exposure>.<file_type>} with spaces replaced by
#' underscores.
#'
#' @return Invisibly returns \code{object} so calls can be chained.
#'
#' @examples
#' \donttest{
#' data("fi_49item")
#' input1   <- harmonize_mr_data(df = fi_49item)
#' outcome1 <- run_mr_analysis(MR_input_data = input1)
#'
#' gwas_plots <- GWAS_forest(MR_input_data = input1, report_form = "Beta")
#' mr_plots   <- MR_forest(summary_df = outcome1, effect = "Beta")
#'
#' # Inspect what is stored; for fi_49item both objects print:
#' #   [1] fi_49item :: Zn
#'
#' # Retrieve the exact outcome/exposure labels stored in the object
#' gwas_out <- gwas_plots@outcomes[1]   # "fi_49item"
#' gwas_exp <- gwas_plots@exposures[1]  # "Zn"
#' mr_out   <- mr_plots@outcomes[1]     # "fi_49item"
#' mr_exp   <- mr_plots@exposures[1]    # "Zn"
#'
#' # Export all instrument-level plots as PDF (commented — writes to disk)
#' export_forest_plots(gwas_plots, save_dir = tempdir())
#'
#' # Export one outcome only
#' export_forest_plots(gwas_plots, save_dir = tempdir(), outcome = gwas_out)
#'
#' # Export one exposure only
#' export_forest_plots(mr_plots, save_dir = tempdir(), exposure = mr_exp)
#'
#' # Export one specific pair
#' export_forest_plots(mr_plots, save_dir = tempdir(), outcome = mr_out, exposure = mr_exp)
#' }
#' @importFrom methods setGeneric setMethod
#' @importFrom ggplot2 ggsave
#' @export
setGeneric("export_forest_plots",
           function(object, save_dir = tempdir(), file_type = "pdf",
                    width = 8, height = 6, dpi = 300,
                    outcome = NULL, exposure = NULL)
             standardGeneric("export_forest_plots"))

# Shared filtering + saving logic used by both methods below
.export_forest_impl <- function(object, save_dir, file_type, width, height,
                                dpi, outcome, exposure, fname_prefix) {

  file_type <- match.arg(file_type, c("pdf", "png", "jpeg", "tiff"))

  # ---- build index of plots to export ----
  keys      <- names(object@plots)
  outcomes  <- object@outcomes
  exposures <- object@exposures

  idx <- seq_along(keys)

  if (!is.null(outcome)) {
    idx <- idx[outcomes[idx] == outcome]
  }
  if (!is.null(exposure)) {
    idx <- idx[exposures[idx] == exposure]
  }

  if (length(idx) == 0) {
    message("No plots matched the supplied outcome/exposure filter.")
    return(invisible(object))
  }

  # ---- save each selected plot ----
  for (i in idx) {
    safe_out <- gsub("\\s+", "_", outcomes[i])
    safe_exp <- gsub("\\s+", "_", exposures[i])
    fname    <- paste0(fname_prefix, safe_out, "_", safe_exp, ".", file_type)
    filepath <- file.path(save_dir, fname)

    ggplot2::ggsave(
      filename = filepath,
      plot     = object@plots[[i]],
      width    = width,
      height   = height,
      dpi      = dpi
    )
    message("Saved: ", filepath)
  }

  invisible(object)
}

#' @describeIn export_forest_plots Method for \code{GWASForestPlots}.
#' @importFrom ggplot2 ggsave
#' @export
setMethod("export_forest_plots", "GWASForestPlots",
          function(object, save_dir = tempdir(), file_type = "pdf",
                   width = 8, height = 6, dpi = 300,
                   outcome = NULL, exposure = NULL) {
            .export_forest_impl(object, save_dir, file_type, width, height,
                                dpi, outcome, exposure,
                                fname_prefix = "Instrument_forest_")
          })

#' @describeIn export_forest_plots Method for \code{MRForestPlots}.
#' @export
setMethod("export_forest_plots", "MRForestPlots",
          function(object, save_dir = tempdir(), file_type = "pdf",
                   width = 8, height = 6, dpi = 300,
                   outcome = NULL, exposure = NULL) {
            .export_forest_impl(object, save_dir, file_type, width, height,
                                dpi, outcome, exposure,
                                fname_prefix = "Method_forest_")
          })
