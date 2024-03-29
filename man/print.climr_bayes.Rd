% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.climr_bayes.R
\name{print.climr_bayes}
\alias{print.climr_bayes}
\title{Print method for an \code{"climr_bayes"} object}
\usage{
\method{print}{climr_bayes}(x, digits = max(3L, getOption("digits") -
  3L), ...)
}
\arguments{
\item{x}{an object of class \code{c("climr_bayes", "climr")}, a result of a
call to \code{\link{fit_bayes}}.}

\item{digits}{The argument \code{digits} to \code{\link{print.default}}.}

\item{...}{Additional arguments.  None are used in this function.}
}
\value{
The argument \code{x}, invisibly, as for all
  \code{\link[base]{print}} methods.
}
\description{
\code{print} method for class \code{c("climr_bayes", "climr")}.
}
\details{
Prints the original call to \code{\link{fit_bayes}} and the
  estimates from Bayesian inference for a two-way random-effects ANOVA model
  using Markov chain Monte Carlo (MCMC) techniques.
}
\section{Examples}{
 See the examples in \code{\link{fit_bayes}}.
}

\seealso{
\code{\link{fit_bayes}} for estimation from Bayesian inference for a
  two-way random-effects ANOVA model using Markov chain Monte Carlo (MCMC)
  techniques.

\code{\link{summary.climr_bayes}}: \code{summary} method for class
  \code{"climr_bayes"}.
}
