#' Print method for objects of class \code{"summary.climr_bayes"}
#'
#' \code{print} method for an object \code{x} of class
#' \code{"summary.climr_bayes"}.
#'
#' @param x An object of class \code{"summary.climr_bayes"}, a result of a call
#'   to \code{\link{summary.climr_bayes}}.
#' @param ... Additional arguments passed on to \code{\link{print.default}}.
#'
#' @details Prints the call and the numeric matrix \code{x$matrix} returned from
#'   \code{\link{summary.climr_bayes}}.
#'
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#'
#' @seealso \code{\link{fit_bayes}} for estimations from Bayesian inference for a
#'   two-way random-effects ANOVA model using Markov chain Monte Carlo (MCMC)
#'   techniques,.
#' @seealso \code{\link{print.climr_bayes}} \code{print} method for class
#'   \code{"climr_bayes"}.
#' @seealso \code{\link{summary.climr_bayes}} \code{summary} method for class
#'   \code{"climr_bayes"}.
#' @section Examples: See the examples in \code{\link{fit_bayes}}.
#'
#' @export
print.summary.climr_bayes <- function(x, ...) {
  if (!inherits(x, "summary.climr_bayes")) {
    stop("use only with \"summary.climr_bayes\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  print(x$matrix, quote = FALSE, ...)
  if (!is.null(x$warning)) {
    cat("\n")
    cat(x$warning)
  }

  invisible(x)
}
