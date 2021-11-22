#' Descriptive plot for Bayesian Analysis
#'
#' \code{\link[graphics]{plot}} method for an object of class
#' \code{"climr_bayes"},following
#' \href{http://dx.doi.org/10.1175/JCLI-D-14-00265.1}{Northrop and Richard
#' (2014)}. Here is two types of plots that are summaries of the posterior
#' distributions for the super-population standard deviations of the four
#' random effects (\code{which = 1, super_pop = TRUE}), summaries of the posterior
#' distributions for the finite-population standard deviations of the four
#' random effects (\code{which = 1, super_pop = FALSE}) and their posterior distrbutions
#' based on the half-Cauchy prior distrbutions (\code{which = 2, super_pop = TRUE}).
#'
#' @param data An object of class \code{"climr_bayes"} form
#'   \code{\link[climr]{fit_bayes}}, containing the \code{brm_object} that is a
#'   result of Bayesian inference (i.e. \code{\link[brms]{brm}}).
#' @param which The choice of types of plots (default to be 1). See also in
#'   Details.
#' @param super_pop The choice of the types of standard deviation. As
#'  \code{super_pop = TRUE}, it shows the super-population standard deviation.
#'  Otherwise, it shows the finite-population standard deviation.
#' @param prob The probability mass to include in the inner interval for
#'   \code{\link[bayesplot]{mcmc_intervals}}, and the default is 0.5 (i.e.50\%
#'   confidence interval)
#' @param prob_outer The probability mass to include in the outer intervalfor
#'   \code{\link[bayesplot]{mcmc_intervals}}, and the default is 0.95 (i.e.95\%
#'   confidence interval)
#' @param point_est	 The point estimate to show. Either "median" (the default),
#'   "mean", or "none".
#' @param A Setting for the parameters of the half-Cauchy prior distributions of
#'   standard deviations of GCM, RCP and their interaction, i.e. the group-level
#'   effects and of the error standard deviation, that is, the standard
#'   deviation over runs for a given GCM-RCP combination.(defaults to be 0.5)
#' @param col_bar The color of the bar in the histogram.(more information in
#'   \code{\link{par}})
#' @param col_line The color of the line (i.e. curve) in the histgram.(more
#'   information in \code{\link{par}})
#' @param legend_pos Position of the legend, nominally,is in \code{"center"}.
#' @param lwd The line width, a positive number, defaulting to 2. (more
#'   information in \code{\link{par}})
#' @param ... Further arguments to be passed to
#'   \code{\link[bayesplot]{mcmc_intervals}}, \code{\link[graphics]{hist}},
#'   \code{\link[graphics]{curve}} and \code{\link[graphics]{legend}}.
#'
#' @details This function is only applicable for the Bayesian result from CMIP5
#'   data whose class is \code{"climr_bayes"} from
#'   \code{\link[climr]{fit_bayes}}, containing the \code{brm_object} that is
#'   the output directly from the \code{\link[brms]{brm}}. Here are two types of
#'   plots in the function.
#'
#'   The first one is the plot from \code{\link[bayesplot]{mcmc_intervals}} as
#'   \code{which = 1}. This plot creates confidence intervals computed from
#'   posterior distribution. For CMIP5, there are four lines interpreting
#'   \code{sd_GCM_Intercept}, \code{sd_RCP_Intercept},
#'   \code{sd_GCM:RCP_Intercept} and \code{sigma}, corresponding to
#'   \eqn{\sigma_{GCM}}{\sigma_GCM}, \eqn{\sigma_{RCP}}{\sigma_RCP},
#'   \eqn{\sigma_{GCM:RCP}}{\sigma_GCM:RCP} and \eqn{\sigma_{run}}{\sigma_run},
#'   respectively. However, for CMIP6, there are four lines interpreting
#'   \code{sd_GCM_Intercept}, \code{sd_SSP_Intercept},
#'   \code{sd_GCM:SSP_Intercept} and \code{sigma}, corresponding to
#'   \eqn{\sigma_{GCM}}{\sigma_GCM}, \eqn{\sigma_{SSP}}{\sigma_SSP},
#'   \eqn{\sigma_{GCM:SSP}}{\sigma_GCM:SSP} and \eqn{\sigma_{run}}{\sigma_run},
#'   respectively. The dot represents the \code{point_est} and the default is
#'   showing the median position. The thick lines are corresponding to the 50\%
#'   confidence interval which could be setted in \code{prob} and the thin lines
#'   are corresponding to the 95\% confidence interval which could also be setted
#'   in \code{prob_outer}. But there is one thing that the value of \code{prob}
#'   should not close to that of \code{prob_outer}. At the same time with
#'   \code{which = 1}, if \code{super_pop = TRUE}, those values are corresponding
#'   to the super-population standard deviations. Otherwise, if \code{super_pop = FALSE},
#'   those values are from the finite-population standard deviations.
#'
#'   If choosing \code{which = 2}, it would show the plot from
#'   \code{\link[graphics]{hist}} creating a histogram of posterior distributions
#'   of \eqn{\sigma_{GCM}}{\sigma_GCM}, \eqn{\sigma_{RCP}}{\sigma_RCP},
#'   \eqn{\sigma_{GCM:RCP}}{\sigma_GCM:RCP} and \eqn{\sigma_{run}}{\sigma_run} for
#'   CMIP5 or of \eqn{\sigma_{GCM}}{\sigma_GCM}, \eqn{\sigma_{SSP}}{\sigma_SSP},
#'   \eqn{\sigma_{GCM:SSP}}{\sigma_GCM:SSP} and \eqn{\sigma_{run}}{\sigma_run}
#'   for CMIP6, both based on the half-Cauchy prior distribution with parameter
#'   \code{A} that is interpreted by the curve from \code{\link[graphics]{curve}}.
#'   From this plot, we could see the fitness of prior distrbuiton with the value
#'   of \code{A} to the posterior distribution. With a suitable value of \code{A}
#'   in the half-Cauchy prior distrbution, there is a hight probability on the
#'   realistic range. But the prior has a "long tail" for preventing an undue
#'   influence from the prior if the data suggest a larger than anticipated value
#'   of \eqn{\sigma}. Besides, the \code{which = 2} is only working for
#'   \code{super_pop = TRUE}. As the \code{which = 2} and \code{super_pop = FALSE}
#'   are choosing at the same time, it would stop.
#'
#'   The relatice mathematics in this
#'   pacakge are seen in \code{vignette("climr", package = "climr")} for an
#'   overview of the package.
#'
#' @return Nothing is returned.
#'
#' @references Northrop, P.J. and R.E. Chandler (2014). Quantifying Sources of
#'   Uncertainty in Projections of Future Climate. \emph{Journal of Climate},
#'   \strong{27}, 8793-8808. \url{https://doi.org/10.1175/JCLI-D-14-00265.1}
#'
#' @seealso \code{\link[climr]{cmip5_temp1}} and
#'   \code{\link[climr]{cmip5_temp2}} for CMIP5 data.
#' @seealso \code{\link[climr]{cmip6_pr}} and \code{\link[climr]{cmip6_tas}} for
#'   CMIP6 data.
#' @seealso \code{\link[climr]{fit_bayes}} for operation of bayesian inference.
#' @seealso \code{\link[bayesplot]{mcmc_intervals}} for bayesplot of the posterior
#' distrbutions
#'
#' @examples
#' x <- fit_bayes(cmip5_temp1)
#' y <- fit_bayes(cmip6_pr)
#' ### Plot method [1] for posterior distribution with super-population standard deviation ###
#' plot(x)
#' plot(y)
#' ### Plot method [1] for posterior distribution with finite-population standard deviation###
#' plot(x, super_pop = FALSE)
#' plot(y, super_pop = FALSE)
#' ### Plot method [2] for prior distribution ###
#' plot(x, which = 2)
#' plot(y, which = 2)
#'
#' @export
plot.climr_bayes <- function(data, which = 1, super_pop = TRUE,  ...) {

  ### Plot method [1] for posterior distribution with super-population standard deviation ###

  plot1 <- function(x, prob = 0.5, ..., prob_outer = 0.95, point_est = c("median")){
    if (inherits(data, "cmip5")) {
      bayesplot::mcmc_intervals(x, prob = prob, prob_outer = prob_outer,
                                point_est = point_est, ...)

    } else if (inherits(data, "cmip6")) {
      bayesplot::mcmc_intervals(x, prob = prob, prob_outer = prob_outer,
                                point_est = point_est, ...)
    }
  }

  ### Plot method [2] for prior distribution ###

  plot2 <- function(y, A, ..., col_bar = "grey", col_line = "red",
                    legend_pos = "topright", lwd = 2){
    my_plot <- function(x, ..., A, xlab, lwd = 2, title) {
      graphics::hist(x, prob = TRUE, main = NULL, xlab = xlab, col = col_bar, ...)
      graphics::curve(2 * dcauchy(x, scale = A), add = TRUE, lwd = lwd, col = col_line, ...)
    }
    # Text for the legends
    a <- paste("half Cauchy (", A, ") prior")
    legend_text <- c(a, "posterior")
    # Produce the plot and add the legend
    # Change mfrow but save the old par() settings in oldpar, so that we can restore them when we exit the function.
    graphics::par(mfrow = c(2, 2))
    my_plot(x = y[,1], A = A, xlab = "standard deviation of GCM", ...)
    # graphics::legend(x = legend_pos, legend = legend_text, pch = c(NA, 15), pt.cex = 2,
    #                 lwd = lwd, col = c(col_line, col_bar), lty = c(1,0), bty='n',
    #                 xpd = TRUE, ...)

    if (inherits(data, "cmip5")) {
      my_plot(x = y[,2], A = A, xlab = "standard deviation of RCP", ...)
      graphics::legend(x = legend_pos, legend = legend_text, pch = c(NA, 15), pt.cex = 2,
                       lwd = lwd, col = c(col_line, col_bar), lty = c(1,0), bty='n',
                       xpd = TRUE, ...)
      my_plot(x = y[,3], A = A, xlab = "standard deviation of GCM:RCP", ...)

    } else if (inherits(data, "cmip6")) {
      my_plot(x = y[,2], A = A, xlab = "standard deviation of SSP",
              legend_text = legend_text, legend_pos = legend_pos, ...)
      graphics::legend(x = legend_pos, legend = legend_text, pch = c(NA, 15), pt.cex = 2,
                       lwd = lwd, col = c(col_line, col_bar), lty = c(1,0), bty='n',
                       xpd = TRUE, ...)
      my_plot(x = y[,3], A = A, xlab = "standard deviation of GCM:SSP", ...)

    }

    my_plot(x = y[,4], A = A, xlab = "standard deviation of run", ...)

  }

  ### Choice of the types of plot ###

  posterior <- as.matrix(data$brm_object)[,c(2,4,3,5)]
  if (inherits(data, "cmip5")) {
    colnames(posterior) <- c("sigma_GCM", "sigma_RCP", "sigma_GCM_RCP", "sigma_run")

  } else if (inherits(data, "cmip6")) {
    colnames(posterior) <- c("sigma_GCM", "sigma_SSP", "sigma_GCM_SSP", "sigma_run")
  }

  finite_sd <- data$finite_sd

  if (inherits(data, "cmip5")) {
    colnames(finite_sd) <- c("s_GCM", "s_RCP", "s_GCM_RCP", "s_run")

  } else if (inherits(data, "cmip6")) {
    colnames(finite_sd) <- c("s_GCM", "s_SSP", "s_GCM_SSP", "s_run")
  }

  if(which == 1 & super_pop == TRUE){ plot1(x = posterior,...) }

  else if(which == 1 & super_pop == FALSE){ plot1(x = finite_sd,...) }

  else if(which == 2 & super_pop == TRUE) { plot2(y = posterior, A = data$brm_object$stanvars$A$sdata, ...) }

  else if(which == 2 & super_pop == FALSE){
    stop("\" which = 2\" only work with  \"super_pop = TRUE\"")
   }
}
