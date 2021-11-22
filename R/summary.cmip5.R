#' Descriptive summary for CMIP5 data
#'
#' \code{summary} method for an object of class \code{"cmip5"},following
#' \href{http://dx.doi.org/10.1175/JCLI-D-14-00265.1}{Northrop and Richard
#' (2014)}
#'
#' @param data An object of class \code{"cmip5"}, a result of a call to
#'   \code{\link{nnt}}.
#' @param ... Further arguments to be passed to \code{\link[base]{summary}}
#'
#' @details This function is only applicable for the CMIP5 data containing the
#'   index, GCM and RCP; Users could obtain an organized table to summary the
#'   number of the runs for each combinations of GCMs and RCPs(i.e. RCP2.6,
#'   RCP4.5, RCP6.0 and RCP8.5).
#'
#' @return return the table \code{k} that contains the overall organised table
#'   for numbers of run for each combination of GCMs and RCPs
#'
#' @references Northrop, P.J. and R.E. Chandler (2014). Quantifying Sources of
#'   Uncertainty in Projections of Future Climate. \emph{Journal of Climate},
#'   \strong{27}, 8793-8808. \url{https://doi.org/10.1175/JCLI-D-14-00265.1}
#'
#' @seealso \code{\link[climr]{cmip5_temp1}} and \code{\link[climr]{cmip5_temp2}}
#'  for CMIP5 data.
#'
#' @examples
#' summary(cmip5_temp1)
#'
#' @export
summary.cmip5 <- function(data,...){
  if (!inherits(x, "cmip5")) {
    stop("use only with \"cmip5\" objects")
  }
  k <- table(data[, c("GCM", "RCP")])
  total <- colSums(k[,1:ncol(k)])
  k <- rbind(k, total)
  row.names(k)[nrow(k)] <- "Total Cases"
  return(k)
}

