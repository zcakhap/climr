#' Descriptive plot for CMIP6 data
#'
#' \code{plot} method for an object of class \code{"cmip6"},following
#' \href{http://dx.doi.org/10.1175/JCLI-D-14-00265.1}{Northrop and Richard
#' (2014)}
#'
#' @param data An object of class \code{"cmip6"}, a result of a call to
#'   \code{\link{nnt}}.
#' @param legend_pos Position of the legend, nominally,is in \code{"topright"}.
#' @param legend_text Text of the legend, nominally, be in the order of
#'   \code{"ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp460",
#'   "ssp534-over", "ssp585"}.
#' @param region Region of the world: 0 for global, 1-22 for regions.(default 0)
#' @param period Time horizon: 1 for mid 21st century (2020-2049); 2 for late
#'   21st century (2070-2099).
#' @param ... Further arguments to be passed to \code{\link[graphics]{plot}} or
#'   \code{\link[graphics]{legend}}.
#'
#' @details This function is only applicable for the CMIP6 data containing the
#'   index, GCM and SSP; Users could adjust the plot whatever they want, like
#'   the \code{pch} and \code{col}. Although there are some arguments to avoid
#'   the systemic misunderstanding between the \code{\link[graphics]{plot}} and
#'   \code{\link[graphics]{legend}}, it is still possible to create warning
#'   messages. But those could be eschewed by careful input.
#'
#'   If users do not have any other further arguments, this function could also
#'   create a satisfied plot for them with classified colors and points.
#'
#'   Colours and types of the points can be set sing an argument \code{col} and
#'   \code{pch}.If a variable is wrapped then the default plotting limits are
#'   set using the corresponding values in \code{data$SSP}.
#'
#' @return Nothing is returned.
#'
#' @references Northrop, P.J. and R.E. Chandler (2014). Quantifying Sources of
#'   Uncertainty in Projections of Future Climate. \emph{Journal of Climate},
#'   \strong{27}, 8793-8808. \url{https://doi.org/10.1175/JCLI-D-14-00265.1}
#'
#' @seealso \code{\link[climr]{cmip6_pr}} and \code{\link[climr]{cmip6_tas}} for
#'   CMIP6 data.
#'
#' @examples
#' plot(cmip6_pr)
#'
#' @export
plot.cmip6 <- function(x, region = 0, period = 1, legend_pos = "topright",
                       legend_text = levels(data$SSP),
                       ...) {
  data <- x[x$region == region & x$period == period,]
  # Extract the index values and set numeric values for GCM and SSP
  index <- data$index
  gcm <- as.numeric(data$GCM)
  ssp <- as.numeric(data$SSP)
  # Create plot and legend functions to avoid problems with ...
  my_plot <- function(x, y, ..., col = 1:8, pch = 1:8, xlab = "GCM number",
                      ylab = "Index of temperature change / deg C", title) {
    # Use rep_len() to ensure that col and pch have length 4
    # Use [rcp] to generate pch and col for each value in index
    graphics::plot(x = x, y = y, col = rep_len(col, 8)[ssp],
                   pch = rep_len(pch, 8)[ssp], xlab = xlab, ylab = ylab, ...)
  }

  my_legend <- function(x, legend, ncol = 2, cex = 0.9, ..., col = 1:8, pch = 1:8, lwd,
                        lty, xlab, ylab, main, axes, frame.plot, xlim, ylim) {
    graphics::legend(x = x, legend = legend, col = col, pch = pch, ncol = ncol, cex = cex, ...)
  }
  # Produce the plot and add the legend
  my_plot(gcm, index, ...)
  my_legend(x = legend_pos, legend = rep_len(legend_text, 8),...)
  return(invisible())
}
