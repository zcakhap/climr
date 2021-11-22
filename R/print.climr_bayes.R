#' Print method for an \code{"climr_bayes"} object
#'
#' \code{print} method for class \code{c("climr_bayes", "climr")}.
#'
#' @param x an object of class \code{c("climr_bayes", "climr")}, a result of a
#'   call to \code{\link{fit_bayes}}.
#' @param digits The argument \code{digits} to \code{\link{print.default}}.
#' @param ... Additional arguments.  None are used in this function.
#'
#' @details Prints the original call to \code{\link{fit_bayes}} and the
#'   estimates from Bayesian inference for a two-way random-effects ANOVA model
#'   using Markov chain Monte Carlo (MCMC) techniques.
#'
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#'
#' @seealso \code{\link{fit_bayes}} for estimation from Bayesian inference for a
#'   two-way random-effects ANOVA model using Markov chain Monte Carlo (MCMC)
#'   techniques.
#' @seealso \code{\link{summary.climr_bayes}}: \code{summary} method for class
#'   \code{"climr_bayes"}.
#' @section Examples: See the examples in \code{\link{fit_bayes}}.
#'
#' @export
print.climr_bayes <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  #checking the class of x
  if (!inherits(x, "climr")) {
    stop("use only with \"climr\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Estimates from the Bayesian inference:\n")

  print.default(format(t(x$coefs), digits = digits), print.gap = 2L,
                quote = FALSE)

  return(invisible(x))
}
