#' Descriptive plot for CMIP6 data
#'
#' \code{summary} method for an object of class \code{"cmip6"},following
#' \href{http://dx.doi.org/10.1175/JCLI-D-14-00265.1}{Northrop and Richard
#' (2014)}
#'
#' @param data An object of class \code{"cmip6"}, a result of a call to
#'   \code{\link{nnt}}.
#' @param region Region of the world: 0 for global, 1-22 for regions.(default 0)
#' @param period Time horizon: 1 for mid 21st century (2020-2049); 2 for late
#'   21st century (2070-2099).
#' @param ... Further arguments to be passed to \code{\link[base]{summary}}
#'
#' @details This function is only applicable for the CMIP6 data containing the
#'   index, GCM and SSP; Users could obtain an organized table to summary the
#'   number of the runs for each combinations of GCMs and SSPs in a specific region
#'   and specific time horizon.
#'
#' @return return the table \code{k} that contains the overall organised table
#'   for numbers of run for each combination of GCMs and SSPs in a specific region
#'   and specific time horizon
#'
#' @references Northrop, P.J. and R.E. Chandler (2014). Quantifying Sources of
#'   Uncertainty in Projections of Future Climate. \emph{Journal of Climate},
#'   \strong{27}, 8793-8808. \url{https://doi.org/10.1175/JCLI-D-14-00265.1}
#'
#' @seealso \code{\link[climr]{cmip6_pr}} and \code{\link[climr]{cmip6_tas}} for
#'   CMIP6 data.
#'
#' @examples
#' summary(cmip6_pr)
#'
#' @export
summary.cmip6 <- function(x, region = 0, period = 1, ...){
  if (!inherits(x, "cmip6")) {
    stop("use only with \"cmip6\" objects")
  }
  data <- x[x$region == region & x$period == period,]
  k <- table(data[, c("GCM", "SSP")])
  total <- colSums(k[,1:ncol(k)])
  k <- rbind(k, total)
  row.names(k)[nrow(k)] <- "Total Cases"
  return(k)
}

