#' Summary method for an \code{"climr_bayes"} object
#'
#' \code{summary} method for class \code{"climr_bayes"}
#'
#' @param x an object of class \code{"climr_bayes"}, a result of a call to
#'   \code{\link{fit_bayes}}.
#' @param digits An integer. Used for number formatting with
#'   \code{\link[base]{summary}}.
#' @param ... Additional arguments.  None are used in this function.
#'
#' @return Returns an object (a list) of class \code{"summary.climr_bayes"}
#'   containing the list element \code{x$call} and a numeric matrix
#'   \code{matrix} giving, for all five parameters, their estimators, the
#'   corresponding estimated standard errors (Std. Error) and 95\% confidence
#'   intervals i.e. the values subtracted from the raw estimate.
#'
#' @seealso \code{\link{fit_bayes}} for estimations from Bayesian inference for
#'   a two-way random-effects ANOVA model using Markov chain Monte Carlo (MCMC)
#'   techniques,.
#' @seealso \code{\link{print.climr_bayes}} \code{print} method for class
#'   \code{"climr_bayes"}.
#' @section Examples: See the examples in \code{\link{fit_bayes}}.
#'
#' @export
summary.climr_bayes <- function(x, digits = max(3, getOption("digits") - 3L),
                               ...) {
  # checking the class of x
  if (!inherits(x, "climr_bayes")) {
    stop("use only with \"climr_bayes\" objects")
  }

  res <- x["call"]
  upper <- signif(x$ci[,1], digits = digits)
  lower <- signif(x$ci[,2], digits = digits)
  new_ci<- paste("(", upper, ",",  lower, ")")

  res$matrix <- cbind(`Estimate` = signif(x$coefs, digits = digits),
                      `Std. Error` = signif(x$se, digits = digits),
                      `95% Confidence Interval` = new_ci)
  colnames(res$matrix) <- c("Estimates","Std.Error","95% Confidence Intervals")
  class(res) <- "summary.climr_bayes"
  return(res)
}
