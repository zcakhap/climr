#' Descriptive plot for CMIP5 data
#'
#' \code{plot} method for an object of class \code{"cmip5"},following
#' \href{http://dx.doi.org/10.1175/JCLI-D-14-00265.1}{Northrop and Richard
#' (2014)}
#'
#' @param data An object of class \code{"cmip5"}, a result of a call to
#'   \code{\link{nnt}}.
#' @param legend_pos Position of the legend, nominally,is in \code{"topright"}.
#' @param legend_text Text of the legend, nominally, be in the order of
#'   \code{"RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"}.
#' @param ... Further arguments to be passed to \code{\link[graphics]{plot}} or
#'   \code{\link[graphics]{legend}}.
#'
#' @details This function is only applicable for the CMIP5 data containing the
#'   index, GCM and RCP; Users could adjust the plot whatever they want, like
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
#'   set using the corresponding values in \code{data$RCP}.
#'
#' @return Nothing is returned.
#'
#' @references Northrop, P.J. and R.E. Chandler (2014). Quantifying
#'   Sources of Uncertainty in Projections of Future Climate.
#'   \emph{Journal of Climate}, \strong{27}, 8793-8808.
#'   \url{https://doi.org/10.1175/JCLI-D-14-00265.1}
#'
#' @seealso \code{\link[climr]{cmip5_temp1}} and \code{\link[climr]{cmip5_temp2}}
#'  for CMIP5 data.
#'
#' @examples
#' plot(cmip5_temp1)
#'
#' @export
plot.cmip5 <- function(data, legend_pos = "topright",
                       legend_text = levels(data$RCP),
                       ...) {
  # Extract the index values and set numeric values for GCM and RCP
  index <- data$index
  gcm <- as.numeric(data$GCM)
  rcp <- as.numeric(data$RCP)
  # Create plot and legend functions to avoid problems with ...
  my_plot <- function(x, y, ..., col = 1:4, pch = 1:4, xlab = "GCM number",
                      ylab = "Index of temperature change / deg C", title) {
    # Use rep_len() to ensure that col and pch have length 4
    # Use [rcp] to generate pch and col for each value in index
    graphics::plot(x = x, y = y, col = rep_len(col, 4)[rcp],
                   pch = rep_len(pch, 4)[rcp], xlab = xlab, ylab = ylab, ...)
  }

  my_legend <- function(x, legend, ..., ncol = 2, cex = 0.9, col = 1:4, pch = 1:4,
                        lwd,lty, xlab, ylab, main, axes, frame.plot, xlim, ylim) {
      graphics::legend(x = x, legend = legend, col = col, pch = pch, cex = cex,
                       ncol = ncol, ...)
  }
  # Produce the plot and add the legend
  my_plot(gcm, index, ...)
  my_legend(x = legend_pos, legend = rep_len(legend_text, 4), ...)
  return(invisible())
}
