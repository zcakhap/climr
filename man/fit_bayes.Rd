% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_bayes.R
\name{fit_bayes}
\alias{fit_bayes}
\title{Two-way Hierarchical Analysis of Variance}
\usage{
fit_bayes(data, region = 0, period = 1, sd_mu = 1e+06, A = 0.5,
  chains = 5, iter = 3000, refresh = 0, adapt_delta = 0.999,
  max_treedepth = 15, seed = 20200705, ...)
}
\arguments{
\item{data}{A dataframe inheriting from class \code{"cmip5"}, containing
variables \code{index}, \code{GCM}, \code{RCP} and \code{run}.  See
\code{\link{cmip5_temp1}}.}

\item{region}{Region of the world specific for CMIP6: 0 for global, 1-22 for
regions.(default 0)}

\item{period}{Time horizon: 1 for mid 21st century (2020-2049); 2 for late
21st century (2070-2099).}

\item{sd_mu}{Setting for the standard deviation of the normal prior
distribution of \eqn{\mu}, i.e population-level effect.}

\item{A}{Setting for the parameters of the half-Cauchy prior distributions of
standard deviations of GCM, RCP and their interaction, i.e. the group-level
effects and of the error standard deviation, that is, the standard
deviation over runs for a given GCM-RCP combination.(defaults to be 0.5)}

\item{chains}{The number of Markov chains, with default value 5.(also see
\code{\link[brms]{brms}}).}

\item{iter}{The number of total iterations per chain (including warmup;
defaults to 3000).(also see \code{\link[brms]{brms}}).}

\item{refresh}{Set \code{refresh = 0} to turn off printing the actual
sampling progress.(also see \code{\link[rstan]{stan}}).}

\item{adapt_delta}{one of parameters in \code{control} that is a named list
of parameters to control the sampler's behavior,including
\code{max_treedepth} as well. Its default value is 0.999.For a
comprehensive overview see \code{\link[rstan]{stan}}.}

\item{max_treedepth}{A maximum depth parameter in \code{control} that is a
named list of parameters to control the sampler's behavior,including
\code{adapt_delta} as well. Its default value is 15. For a comprehensive
overview see \code{\link[rstan]{stan}}.}

\item{seed}{The seed for random number generation to make results
reproducible (default to be 20200705 setted by randomly).}

\item{...}{Further arguments to be passed to \code{\link[brms]{brm}}.}
}
\value{
An object (a list) named \code{res} containing the following
  components. \item{call}{the call to \code{fit_bayes}} \item{brm_object
  }{the result of fitted Bayesian model using the MCMC method from
  \code{\link[brms]{brm}}} \item{coefs}{containing the estimated coefficients
  for \eqn{\mu}, i.e. population-effect factor, and the estimated super-population
  standard deviation of GCM, RCP, GCM:RCP and run for CMIP5 or of GCM, SSP, GCM:SSP
  and run for CMIP6, i.e. group-effect factors} \item{se}{containing the
  standard errors for \eqn{\mu}, i.e. population-effect factor, and GCM, RCP,
  GCM:RCP and run for CMIP5 or of GCM, SSP, GCM:SSP and run for CMIP6, i.e.
  group-effect factors} \item{coefs}{containing the 95\% confidence intervals
  for \eqn{\mu}, i.e. population-effect factor, and GCM, RCP, GCM:RCP and run
  for CMIP5 or of GCM, SSP, GCM:SSP and run for CMIP6, i.e. group-effect
  factors}\item{finite_sd}{containing the finite-population standard deviations of
  GCM, RCP, GCM:RCP and run for CMIP5 or of GCM, SSP, GCM:SSP and run for CMIP6,
  i.e. group-effect factors}
}
\description{
Bayesian inference for a two-way random-effects ANOVA model using Markov
chain Monte Carlo (MCMC) techniques, following
\href{http://dx.doi.org/10.1175/JCLI-D-14-00265.1}{Northrop and Chandler
(2014)}.
}
\details{
Those mathematics in this pacakge are seen in
  \code{vignette("climr", package = "climr")} and
  \code{vignette("TwHAoV_bayes", package = "climr")} for an overview of the
  package.

  The data would fit in a Bayesian generalized multivariate multilevel model
  containing three factors which are GCM, RCP and the interaction between GCM
  and RCP for CMIP 5 or which are GCM, SSP and the interaction between GCM
  and SSP for CMIP6. Their prior distributions of standard deviations are
  considered to be half-Cauchy with parameter \code{A} that defaults to be
  0.5. This is same as the prior distribution of residual, i.e. run. For the
  intercept, its standard deviation has normal prior distribution with zero
  mean and \code{sd_mu} standard deviation. \code{sd_mu} defaults to be 1e6.

  After setting the prior distributions, we use MCMC techniques to establish
  a Bayesian model. In addition to choosing the number of chains and the
  number of total iteration per chain, users can control the behavior of the
  NUTS sampler, by using the \code{control} argument. The most important
  reason to use control is to decrease (or eliminate at best) the number of
  divergent transitions that cause a bias in the obtained posterior samples.
  Whenever you see the warning "There were x divergent transitions after
  warmup." you should really think about increasing \code{adapt_delta}.
  Increasing \code{adapt_delta} will slow down the sampler but will decrease
  the number of divergent transitions threatening the validity of your
  posterior samples. Thus, we normally choose the \code{adapt_delta} to be
  0.999. In addition, there would be some warning messages to show that there
  were a lot transitions after warmup that exceeded the maximum treedepth.By
  increasing the parameter of \code{max_treedepth} with default value 15 in
  \code{control}, this problem could solved. Besides, if the Bulk Effective
  Sample Size(ESS) is too low, indicating posterior means and medians may be
  unreliable, we could increase the number of iteration to solve the problem.
  In general, the higher the ESS the better.We recommend that the bulk-ESS is
  greater than 100 times the number of chains. For example, when running five
  chains, this corresponds to having a rank-normalized effective sample size
  of at least 500. The default value of iteration is 2000 in
  \code{\link[brms]{brm}} but for this function, 2000 is too low so we change
  the default value to 3000. More information see \code{\link[rstan]{stan}}.

  Since the whole sampling process is random, it would cause different
  warning messages. Thus, we should choose a number by \code{seed} to choose
  a specific number before running the function. And according to previous
  description matching with the warning messages, we should do some
  adjustments. Inside the function, we are randomly choosing 20200705 as the
  default number of seed.
}
\examples{
library(brms)
## For cmip5_temp1 ##
x <- fit_bayes(cmip5_temp1)

## For cmip5_temp2 ##
y <- fit_bayes(cmip5_temp2)

#' ## For cmip6_pr ##
a <- fit_bayes(cmip6_pr)
b <- fit_bayes(cmip6_pr, period = 2)

## For cmip5_temp2 ##
c <- fit_bayes(cmip6_tas)
d <- fit_bayes(cmip6_tas, period = 2)

}
\references{
Northrop, P.J. and R.E. Chandler (2014). Quantifying Sources of
  Uncertainty in Projections of Future Climate. \emph{Journal of Climate},
  \strong{27}, 8793-8808. \url{https://doi.org/10.1175/JCLI-D-14-00265.1}
}
\seealso{
\code{\link[climr]{cmip5_temp1}} and
  \code{\link[climr]{cmip5_temp2}} for CMIP5 data.

\code{\link[climr]{cmip6_pr}} and \code{\link[climr]{cmip6_tas}} for
  CMIP6 data.

\code{\link[brms]{brm}} for further arguments.

\code{\link[rstan]{stan}} for \code{control}.

\code{\link[climr]{print.climr_bayes}} for \code{print}.

\code{\link[climr]{summary.climr_bayes}} for \code{summary}.

\code{\link[climr]{print.summary.climr_bayes}} for
  \code{print.summary}.
}
