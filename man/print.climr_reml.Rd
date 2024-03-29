% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.climr_reml.R
\name{print.climr_reml}
\alias{print.climr_reml}
\title{Print method for an \code{"climr_reml"} object}
\usage{
\method{print}{climr_reml}(x, digits = max(3L, getOption("digits") - 3L),
  ...)
}
\arguments{
\item{x}{an object of class \code{c("climr_reml", "climr")}, a result of a
call to \code{\link{fit_reml}}.}

\item{digits}{The argument \code{digits} to \code{\link{print.default}}.}

\item{...}{Additional arguments.  None are used in this function.}
}
\value{
The argument \code{x}, invisibly, as for all
  \code{\link[base]{print}} methods.
}
\description{
\code{print} method for class \code{c("climr_reml", "climr")}.
}
\details{
Prints the original call to \code{\link{fit_reml}} and the estimates
  of a fitted ANOVA model using the restricted maximum likelihood(REML).
}
\section{Examples}{
 See the examples in \code{\link{fit_reml}}.
}

\seealso{
\code{\link{fit_reml}} for estimation of a fitted ANOVA model using
  the restricted maximum likelihood(REML).

\code{\link{summary.climr_reml}}: \code{summary} method for class
  \code{"climr_reml"}.
}
