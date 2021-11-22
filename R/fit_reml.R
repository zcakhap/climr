#' Two-way Hierarchical Analysis of Variance
#'
#' Description fitting a random-effects ANOVA model and finding the standard
#' deviations of estimators by using the restricted maximum likelihood (REML)
#' estimation, following
#' \href{http://dx.doi.org/10.1175/JCLI-D-14-00265.1}{Northrop and Chandler
#' (2014)}
#'
#' @param data A dataframe inheriting from class \code{"cmpi5"}, containing
#'   variables \code{index}, \code{GCM}, \code{RCP} and \code{run}.  See
#'   \code{\link{cmip5_temp1}}.
#' @param control The argument \code{control} to
#'   \code{\link[lme4]{lmerControl}}, normally, using the optimizer named
#'   \code{"Nelder_Mead"}.
#' @param ... Further arguments to be passed to \code{\link[lme4]{lmer}}.
#'
#' @details The data would fit in a two-way random-effects ANOVA model
#'   containing three factors which are GCM, RCP and the interaction between GCM
#'   and RCP. Those factors are all considered to be normally distributed random
#'   variables with mean zero and different standard deviation, like
#'   \eqn{\sigma_{GCM)},\sigma_{RCP},\sigma_{GCM:RCP}}{\sigma_(GCM),\sigma_(RCP),\sigma_(GCM:RCP)}.
#'   In addition, the intercept \eqn{\mu} of this model interprets the overall
#'   mean change in the index and the error term \eqn{\epsilon} are independent
#'   identically distributed random variables with mean zero and standard
#'   deviation \eqn{\sigma_{run}}{\sigma_(run)}.
#'
#'   This function uses the restricted maximum likelihood (REML) method for
#'   estimating those estimators and their corresponding standard errors.
#'
#' @return An object (a list) of class \code{climr_reml} and \code{climr}
#'   containing the following components.
#'
#'   \item{call}{the call of \code{fit_reml}} \item{lmer_object}{the result of
#'   fitted model using the REML method from \code{\link[lme4]{lmer}}}
#'   \item{ests}{containing the estimated coefficients for \eqn{\mu} and the
#'   estimated standard deviations of random-effects factors} \item{ses}{the
#'   standard errors of estimators for fixed-effects and that of estimated
#'   standard deviations of random-effects factors} \item{ci}{containing the
#'   upper limits (i.e 2.5\%) and lower limits (i.e 97.5\% ) of the 95\% confidence
#'   intervals of estimators for fixed-effects and that of estimated standard
#'   deviations of random-effects factors}
#'
#' @references Northrop, P.J. and R.E. Chandler (2014). Quantifying
#'   Sources of Uncertainty in Projections of Future Climate.
#'   \emph{Journal of Climate}, \strong{27}, 8793-8808.
#'   \url{https://doi.org/10.1175/JCLI-D-14-00265.1}
#'
#' @seealso \code{\link[climr]{cmip5_temp1}} and \code{\link[climr]{cmip5_temp2}}
#'  for CMIP5 data.
#' @seealso \code{\link[lme4]{lmer}} for further arguments.
#' @seealso \code{\link[climr]{print.climr_reml}} for \code{print}.
#' @seealso \code{\link[climr]{summary.climr_reml}} for \code{summary}.
#' @seealso \code{\link[climr]{print.summary.climr_reml}} for \code{print.summary}.
#' @examples
#' library(lme4)
#' res <- fit_reml(cmip5_temp1)
#' @export
fit_reml <- function(data, control = lme4::lmerControl(optimizer = "Nelder_Mead"),
                     ...) {
  Call <- match.call()
  # Fit model using REML
  x <- lme4::lmer(index ~ (1 | RCP) + (1 | GCM) + (1 | RCP:GCM),
                  data = data, control = control, ...)
  # Create a list containing the results that we need, and the object
  # returned from lme4::lmer()
  res <- list()
  res$call <- Call
  res$lmer_object <- x
  #
  # Estimates --------
  #
  # Estimate of the fixed-effects (intercept)
  mu <- as.numeric(lme4::fixef(x))
  # Estimates of the random-effects standard deviations
  random_effects <- as.numeric(unlist(arm::sigma.hat(x))[1:4])
  # lmer returns these in a strange order!
  sigma_alpha <- random_effects[3]
  sigma_beta <- random_effects[4]
  sigma_gamma <- random_effects[2]
  sigma_epsilon <- random_effects[1]
  res$ests <- c(mu, sigma_alpha, sigma_beta, sigma_gamma, sigma_epsilon)
  names(res$ests) <- c("mu", "sigma_GCM", "sigma_RCP", "sigma_GCM_RCP",
                       "sigma_run")
  #
  # Estimated standard errors --------
  #
  se_mu <- arm::se.fixef(x)

  x.ML <- lme4:::devfun2(x, useSc=TRUE, signames=FALSE)
  se <- as.data.frame(lme4::VarCorr(x))  ## need ML estimates!
  hh1 <- numDeriv::hessian(x.ML,se[,"sdcor"])
  vv2 <- 2*solve(hh1)  ## 2* converts from log-likelihood to deviance scale
  ses <- sqrt(diag(vv2))  ## get standard errors

  res$ses <- c(se_mu, ses)
  names(res$ses) <- c("mu", "sigma_GCM", "sigma_RCP", "sigma_GCM_RCP",
                       "sigma_run")

  # 95% confidence interval
  ci <- lme4::confint.merMod(x, quiet = TRUE)
  res$ci <- ci[c(5,1,2,3,4),]
  names(res$ci) <- c("mu", "sigma_GCM", "sigma_RCP", "sigma_GCM_RCP",
                      "sigma_run")
  class(res) <- c("climr_reml", "cmip5", "climr")
  return(res)
}
