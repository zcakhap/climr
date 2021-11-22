#' Print method for objects of class \code{"summary.climr_reml"}
#'
#' \code{print} method for an object \code{x} of class \code{"summary.climr_reml"}.
#'
#' @param x An object of class \code{"summary.climr_reml"}, a result of a call to
#'   \code{\link{summary.climr_reml}}.
#' @param ... Additional arguments passed on to \code{\link{print.default}}.
#'
#' @details Prints the call and the numeric matrix \code{x$matrix} returned from
#'   \code{\link{summary.climr_reml}}.
#'
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#'
#' @seealso \code{\link{fit_reml}} for estimation of a fitted ANOVA model using
#'   the restricted maximum likelihood(REML).
#' @seealso \code{\link{print.climr_reml}} \code{print} method for class
#'   \code{"climr_reml"}.
#' @seealso \code{\link{summary.climr_reml}}: \code{summary} method for class
#'   \code{"climr_reml"}.
#' @section Examples: See the examples in \code{\link{fit_reml}}.
#'
#' @export
print.summary.climr_reml <- function(x, ...) {
  if (!inherits(x, "summary.climr_reml")) {
    stop("use only with \"summary.climr_reml\" objects")
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
