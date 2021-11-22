#' Print method for an \code{"climr_reml"} object
#'
#' \code{print} method for class \code{c("climr_reml", "climr")}.
#'
#' @param x an object of class \code{c("climr_reml", "climr")}, a result of a
#'   call to \code{\link{fit_reml}}.
#' @param digits The argument \code{digits} to \code{\link{print.default}}.
#' @param ... Additional arguments.  None are used in this function.
#'
#' @details Prints the original call to \code{\link{fit_reml}} and the estimates
#'   of a fitted ANOVA model using the restricted maximum likelihood(REML).
#'
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#'
#' @seealso \code{\link{fit_reml}} for estimation of a fitted ANOVA model using
#'   the restricted maximum likelihood(REML).
#' @seealso \code{\link{summary.climr_reml}}: \code{summary} method for class
#'   \code{"climr_reml"}.
#' @section Examples: See the examples in \code{\link{fit_reml}}.
#' @export
print.climr_reml <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  #checking the class of x
  if (!inherits(x, "climr")) {
    stop("use only with \"climr\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Estimates from the restricted maximum likelihood(REML):\n")
  coefs <- x$ests
  print.default(format(coefs, digits = digits), print.gap = 2L,
                quote = FALSE)

  return(invisible(x))
}
