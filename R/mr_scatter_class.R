# ==============================================================================
# S4 Class: MRScatterPlots
#
# A typed container for the named list of plot-parameter lists produced by
# plot_mr_scatter(). Each element holds all data and settings needed to
# render one scatter plot. Actual drawing is performed by .draw_scatter_plot()
# at export time, avoiding recordPlot()/replayPlot() device dependencies.
# ==============================================================================

# ------------------------------------------------------------------------------
# Class definition
# ------------------------------------------------------------------------------

#' S4 class to store MR scatter plots
#'
#' An \code{MRScatterPlots} object is returned by \code{plot_mr_scatter()}.
#' It holds one plot-parameter list per outcome-exposure pair, together
#' with the metadata needed to filter and export subsets of plots.
#'
#' @slot plots    Named list of plot-parameter lists. Names follow the
#'   pattern \code{"Outcome::Exposure"}.
#' @slot outcomes Character vector of outcome labels, parallel to \code{plots}.
#' @slot exposures Character vector of exposure labels, parallel to
#'   \code{plots}.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{export_scatter_plots(object, save_dir, file_type, width,
#'     height, outcome, exposure)}}{
#'     Saves plots to disk. Optionally filter by \code{outcome},
#'     \code{exposure}, or both to export a subset.}
#' }
#' @importFrom methods setClass setValidity new
#' @export
setClass("MRScatterPlots", representation(
  plots     = "list",
  outcomes  = "character",
  exposures = "character"
))

setValidity("MRScatterPlots", function(object) {
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
# Internal constructor (not exported — called only by plot_mr_scatter)
# ------------------------------------------------------------------------------

#' @noRd
new_mr_scatter_plots <- function(plots, outcomes, exposures) {
  new("MRScatterPlots",
      plots     = plots,
      outcomes  = outcomes,
      exposures = exposures)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# export generic and method
# ------------------------------------------------------------------------------

#' Export MR Scatter Plots to Disk
#'
#' Saves plots stored in an \code{MRScatterPlots} object to a directory.
#' You can export all plots or filter to a specific outcome, a specific
#' exposure, or a specific outcome-exposure pair.
#'
#' @param object    An \code{MRScatterPlots} object returned by
#'   \code{plot_mr_scatter()}.
#' @param save_dir  Character string; directory to write files into.
#'   Must exist. Defaults to \code{tempdir()}.
#' @param file_type Character string; output format passed to the
#'   corresponding \code{grDevices} device. One of \code{"pdf"},
#'   \code{"png"}, \code{"jpeg"}, or \code{"tiff"}.
#'   Default is \code{"png"}.
#' @param width     Numeric; plot width in inches. Default is \code{8}.
#' @param height    Numeric; plot height in inches. Default is \code{6}.
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
#' File names follow the pattern
#' \code{Scatter_<Outcome>_<Exposure>.<file_type>} with spaces replaced
#' by underscores.
#'
#' @return Invisibly returns \code{object} so calls can be chained.
#'
#' @examples
#' \donttest{
#' data("fi_49item")
#' input1   <- harmonize_mr_data(df = fi_49item)$input_df
#' outcome1 <- run_mr_analysis(MR_input_data = input1)
#' plots    <- plot_mr_scatter(MR_input_data = input1, summary_df = outcome1)
#'
#' # Inspect what is stored; for fi_49item this prints:
#' #   [1] fi_49item :: Zn
#'#' # Retrieve the exact outcome/exposure labels stored in the object
#' out_name <- plots@outcomes[1]   # "fi_49item"
#' exp_name <- plots@exposures[1]  # "Zn"
#'
#' # Export all plots to tempdir() (commented — writes to disk)
#' export_scatter_plots(plots, save_dir = tempdir())
#'
#' # Export one outcome only
#' export_scatter_plots(plots, save_dir = tempdir(), outcome = out_name)
#'
#' # Export one exposure only
#' export_scatter_plots(plots, save_dir = tempdir(), exposure = exp_name)
#'
#' # Export one specific pair
#' export_scatter_plots(plots, save_dir = tempdir(), outcome = out_name, exposure = exp_name)
#' }
#' @importFrom methods setGeneric setMethod
#' @importFrom grDevices pdf png jpeg tiff dev.off
#' @export
setGeneric("export_scatter_plots",
           function(object, save_dir = tempdir(), file_type = "png",
                    width = 8, height = 6,
                    outcome = NULL, exposure = NULL)
             standardGeneric("export_scatter_plots"))

#' @describeIn export_scatter_plots Method for \code{MRScatterPlots}.
#' @importFrom grDevices pdf png jpeg tiff dev.off
#' @export
setMethod("export_scatter_plots", "MRScatterPlots",
          function(object, save_dir = tempdir(), file_type = "png",
                   width = 8, height = 6,
                   outcome = NULL, exposure = NULL) {

            # ---- validate file_type ----
            file_type <- match.arg(file_type, c("pdf", "png", "jpeg", "tiff"))

            # ---- validate save_dir ----
            if (!dir.exists(save_dir))
              stop("'save_dir' does not exist: ", save_dir, call. = FALSE)

            # ---- determine which indices to export based on filter arguments ----
            idx <- seq_along(object@plots)

            if (!is.null(outcome) && !is.null(exposure)) {
              # Case 3: specific pair
              idx <- which(object@outcomes == outcome & object@exposures == exposure)
              if (length(idx) == 0)
                stop("No plot found for outcome = '", outcome,
                     "' and exposure = '", exposure, "'.", call. = FALSE)

            } else if (!is.null(outcome)) {
              # Case 1: one outcome, all its exposures
              idx <- which(object@outcomes == outcome)
              if (length(idx) == 0)
                stop("No plots found for outcome = '", outcome, "'.", call. = FALSE)

            } else if (!is.null(exposure)) {
              # Case 2: one exposure, all its outcomes
              idx <- which(object@exposures == exposure)
              if (length(idx) == 0)
                stop("No plots found for exposure = '", exposure, "'.", call. = FALSE)
            }
            # else: idx stays as all indices

            # ---- open device helper ----
            open_device <- function(filepath) {
              switch(file_type,
                     pdf  = grDevices::pdf( filepath, width = width, height = height),
                     png  = grDevices::png( filepath, width = width, height = height,
                                            units = "in", res = 300),
                     jpeg = grDevices::jpeg(filepath, width = width, height = height,
                                            units = "in", res = 300),
                     tiff = grDevices::tiff(filepath, width = width, height = height,
                                            units = "in", res = 300)
              )
            }

            # ---- export each selected plot ----
            for (i in idx) {
              fname <- paste0(
                "Scatter_",
                gsub("\\s+", "_", object@outcomes[i]),  "_",
                gsub("\\s+", "_", object@exposures[i]), ".", file_type
              )
              filepath <- file.path(save_dir, fname)
              open_device(filepath)
              tryCatch({
                .draw_scatter_plot(object@plots[[i]])
                grDevices::dev.off()
              }, error = function(e) {
                grDevices::dev.off()
                stop("Failed to draw plot for '", object@outcomes[i],
                     " vs ", object@exposures[i], "': ", conditionMessage(e),
                     call. = FALSE)
              })
              message("Saved: ", filepath)
            }

            invisible(object)
          })

# ==============================================================================
# showplot() generic and methods
# ==============================================================================

#' Display MR Plots on Screen
#'
#' Renders plots stored in an \code{MRScatterPlots}, \code{GWASForestPlots},
#' or \code{MRForestPlots} object to the active graphics device. You can
#' display all plots or filter to a specific outcome, exposure, or
#' outcome-exposure pair.
#'
#' @param object An \code{MRScatterPlots}, \code{GWASForestPlots}, or
#'   \code{MRForestPlots} object.
#' @param outcome Optional character string; if supplied, only plots whose
#'   outcome matches this value are displayed.
#' @param exposure Optional character string; if supplied, only plots whose
#'   exposure matches this value are displayed.
#'
#' @details Filtering logic:
#' \itemize{
#'   \item \strong{Both \code{outcome} and \code{exposure} supplied} — displays
#'     the single plot for that exact pair.
#'   \item \strong{\code{outcome} only} — displays all exposures for that
#'     outcome.
#'   \item \strong{\code{exposure} only} — displays all outcomes for that
#'     exposure.
#'   \item \strong{Neither supplied} — displays every plot in the object.
#' }
#'
#' @return Invisibly returns \code{object}.
#'
#' @examples
#' \donttest{
#' data("fi_49item")
#' input1   <- harmonize_mr_data(df = fi_49item)$input_df
#' outcome1 <- run_mr_analysis(MR_input_data = input1)
#'
#' # Scatter plots
#' scatter <- plot_mr_scatter(MR_input_data = input1, summary_df = outcome1)
#' showplot(scatter)
#' showplot(scatter, outcome = "fi_49item")
#' showplot(scatter, outcome = "fi_49item", exposure = "Zn")
#'
#' # GWAS forest plots
#' gwas_plots <- GWAS_forest(MR_input_data = input1, report_form = "Beta")
#' showplot(gwas_plots)
#' showplot(gwas_plots, exposure = "Zn")
#'
#' # MR forest plots
#' mr_plots <- MR_forest(summary_df = outcome1, effect = "Beta")
#' showplot(mr_plots)
#' showplot(mr_plots, outcome = "fi_49item", exposure = "Zn")
#' }
#' @importFrom methods setGeneric setMethod
#' @export
setGeneric("showplot",
           function(object, outcome = NULL, exposure = NULL)
             standardGeneric("showplot"))

#' @describeIn showplot Display scatter plots from an \code{MRScatterPlots} object.
#' @export
setMethod("showplot", "MRScatterPlots", function(object, outcome = NULL, exposure = NULL) {
  idx <- .showplot_filter(object, outcome, exposure)
  for (i in idx) {
    .draw_scatter_plot(object@plots[[i]])
  }
  invisible(object)
})

#' @describeIn showplot Display forest plots from a \code{GWASForestPlots} object.
#' @importFrom graphics plot
#' @export
setMethod("showplot", "GWASForestPlots", function(object, outcome = NULL, exposure = NULL) {
  idx <- .showplot_filter(object, outcome, exposure)
  for (i in idx) {
    print(object@plots[[i]])
  }
  invisible(object)
})

#' @describeIn showplot Display forest plots from an \code{MRForestPlots} object.
#' @export
setMethod("showplot", "MRForestPlots", function(object, outcome = NULL, exposure = NULL) {
  idx <- .showplot_filter(object, outcome, exposure)
  for (i in idx) {
    print(object@plots[[i]])
  }
  invisible(object)
})

# Internal helper: resolve filtering to a vector of indices
#' @noRd
.showplot_filter <- function(object, outcome, exposure) {
  idx <- seq_along(object@plots)
  if (!is.null(outcome) && !is.null(exposure)) {
    idx <- which(object@outcomes == outcome & object@exposures == exposure)
    if (length(idx) == 0)
      stop("No plot found for outcome = '", outcome,
           "' and exposure = '", exposure, "'.", call. = FALSE)
  } else if (!is.null(outcome)) {
    idx <- which(object@outcomes == outcome)
    if (length(idx) == 0)
      stop("No plots found for outcome = '", outcome, "'.", call. = FALSE)
  } else if (!is.null(exposure)) {
    idx <- which(object@exposures == exposure)
    if (length(idx) == 0)
      stop("No plots found for exposure = '", exposure, "'.", call. = FALSE)
  }
  idx
}
