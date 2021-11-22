#' Summary method for an \code{"climr_reml"} object
#'
#' \code{summary} method for class \code{"climr_reml"}
#'
#' @param x an object of class \code{"climr_reml"}, a result of a call to
#'   \code{\link{fit_reml}}.
#' @param digits An integer. Used for number formatting with
#'   \code{\link[base]{summary}}.
#' @param ... Additional arguments.  None are used in this function.
#'
#' @return Returns an object (a list) of class \code{"summary.climr_reml"}
#'   containing the list element \code{x$call} and a numeric matrix
#'   \code{matrix} giving, for all five parameters, their estimators and the
#'   corresponding estimated standard errors (Std. Error) 95\% confidence
#'   intervals i.e. the values subtracted from the raw estimate.
#'
#' @seealso \code{\link{fit_reml}} for estimation of a fitted ANOVA model using
#'   the restricted maximum likelihood(REML).
#' @seealso \code{\link{print.climr_reml}} \code{print} method for class
#'   \code{"climr_reml"}.
#' @section Examples: See the examples in \code{\link{fit_reml}}.
#'
#' @export
summary.climr_reml <- function(x, digits = max(3, getOption("digits") - 3L),
                        ...) {
  # checking the class of x
  if (!inherits(x, "climr_reml")) {
    stop("use only with \"climr_reml\" objects")
  }

  res <- x["call"]
  upper <- signif(x$ci[,1], digits = digits)
  lower <- signif(x$ci[,2], digits = digits)
  new_ci<- paste("(", upper, ",",  lower, ")")

  res$matrix <- cbind(`Estimate` = signif(x$ests, digits = digits),
                      `Std. Error` = signif(x$ses, digits = digits),
                      `95% Confidence Interval` = new_ci)
  colnames(res$matrix) <- c("Estimates","Std.Error","95% Confidence Intervals")

  class(res) <- "summary.climr_reml"
  return(res)
}
