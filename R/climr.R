#' climr: Analysis of Climate Variability
#'
#' @description The R package named \code{climr} is serving for analysis of
#'   variability of climate. There are four sets of data in the package,
#'   \code{cmip5_temp1}, \code{cmip5_temp2}, \code{cmip6_tas} and \code{cmip6_pr},
#'   respectively. By analysing the CMIP5 and CMIP6 climate data, we uses the
#'   two-way random-effected ANOVA model for analysing the variability of climate
#'   and finding the main source of climate uncertainty, following the struture of
#'   \href{http://dx.doi.org/10.1175/JCLI-D-14-00265.1}{Northrop and Chandler
#'   (2014)} .This could lead to a great understanding of past, present and
#'   future climate change and variability in a multi-model framework.
#'
#'   We use the two-way random-effected ANOVA model for analysing the CMIP5 and
#'   CMIP6 to find out the climate variability. For CMIP5, we are looking for the
#'   estimates of \eqn{\mu},\eqn{\sigma_{GCM}}{\sigma_GCM},
#'   \eqn{\sigma_{RCP}}{\sigma_RCP}, \eqn{\sigma_{GCM:RCP}}{\sigma_GCM:RCP} and
#'   \eqn{\sigma_{run}}{\sigma_run}. And for  CMIP6 climate, we are looking for the
#'   estimates of \eqn{\mu},\eqn{\sigma_{GCM}}{\sigma_GCM},
#'   \eqn{\sigma_{SSP}}{\sigma_SSP}, \eqn{\sigma_{GCM:SSP}}{\sigma_GCM:SSP} and
#'   \eqn{\sigma_{run}}{\sigma_run}.Those could lead to a great understanding of
#'   past, present and future climate change and variability in a multi-model
#'   framework.
#'
#' @details There are two types of data, CMIP5 and CMIP6, namely. For the CMIP5, the
#'   \code{cmip5_temp1} is showing the indices of global temperature change from
#'   late 20th century (1970-1999) to mid 21st century (2020-2049) and
#'   \code{cmip5_temp2} is showing the indices of global temperature change from
#'   late 20th century (1970-1999) to late 21st century (2069-2098). Both of them
#'   have 270 rows and 4 columns.The first column is the index of temperature change
#'   from one of 38 different General Circulation Models (GCMs) (i.e. in
#'   the second column) under a particular Representative Concentration Pathway
#'   (RCP) (i.e. in the third column) with the number of runs (i.e. in
#'   the forth column). However, for the CMIP6, the \code{cmip6_pr} contains the
#'   indices of percentage precipitation change from late 20th century (1970-1999)
#'   to mid 21st century (2020-2049) and to late 21st century (2070-2099) and the
#'   climate variable is the percipirarion flux(pr). And the \code{cmip6_tas} contains
#'   the indices of temperature change from late 20th century (1970-1999) to mid 21st
#'   century (2020-2049) and to late 21st century (2070-2099) and the climate variable
#'   is surface air temperature(tas). Comparing with the CMIP5, RCP is replaced by
#'   Shared Socio-Economic Pathway (SSP) in CMIP6. RCP has four levels which are
#'   rcp26, rcp45, rcp60 andrcp85 but SSP has 8 levels which are ssp119, ssp126,
#'   ssp245, ssp370, ssp434, ssp460, ssp534-over and ssp585. Besides, for CMIP6, there
#'   are two extra columns (i.e. column 5 and column 6 ) for interpreting the region
#'   and the period.
#'
#'   Besides, we are creating S3 method for demonstrating the data. For CMIP5, we have
#'   \code{\link[climr]{plot.cmip5}} for ploting the graphies and
#'   \code{\link[climr]{summary.cmip5}} for summaring the data into an organised table
#'   from \code{cmip5_temp1} and \code{cmip5_temp2}. Meanwhile, for CMIP6, we have
#'   \code{\link[climr]{plot.cmip6}} for ploting the graphies and
#'   \code{\link[climr]{summary.cmip6}} for summaring the data into an organised table
#'   from \code{cmip6_pr} and \code{cmip6_tas}.
#'
#'   Then we are mainly using the two-way random-effected ANOVA model to analyse
#'   the variability of climate.For example, for CMIP5, we are treated the index of
#'   temperature change to be the \eqn{Y_{ijk}}{Y_(ijk)} corresponding to the \eqn{i}th
#'    GCM, \eqn{j}th RCP and \eqn{k}th run. (e.g. for the CMIP5 data, \eqn{i} varies
#'   from 1 to 38, i.e. the number of GCM; \eqn{j} varies from 1 to 4, i.e. the
#'   number of RCP; and \eqn{k} varies from 0 to 13, i.e. the number of run, to
#'   corresponding the specific GCM and RCP). Besidse, there are one coefficient
#'   (\eqn{\mu}) and five random variables that are \eqn{\alpha_i} for GCM,
#'   \eqn{\beta_j} for RCP, \eqn{\gamma_{ij}}{\gamma_(ij)} for the interaction
#'   between the GCM and RCP and \eqn{\epsilon_{ijk}}{\epsilon_(ijk)} for the run.
#'   Those random varibles are all normally distributed with 0 mean and specific
#'   variance (i.e. \eqn{\sigma_{GCM}}{\sigma_GCM}, \eqn{\sigma_{RCP}}{\sigma_RCP},
#'   \eqn{\sigma_{GCM:RCP}}{\sigma_GCM:RCP} and \eqn{\sigma_{run}}{\sigma_run}).
#'   Thus, we are mainly estimating the \eqn{\mu}, \eqn{\sigma_{GCM}}{\sigma_GCM},
#'   \eqn{\sigma_{RCP}}{\sigma_RCP}, \eqn{\sigma_{GCM:RCP}}{\sigma_GCM:RCP} and
#'   \eqn{\sigma_{run}}{\sigma_run} by those corresponding estimates,
#'   standard errors and 95\% confidence interval. Further, for CMIP6, we are using the
#'   same method but using SSP instead of RCP. Those mathematics in this
#'   pacakge are seen in \code{vignette("climr", package = "climr")} for an
#'   overview of the package.
#'
#'   Within the package, there are two functions representing two difference
#'   inferences to obtain those estimates, following
#'   \href{http://dx.doi.org/10.1175/JCLI-D-14-00265.1}{Northrop and Chandler
#'   (2014)}. One is \code{\link[climr]{fit_reml}} by using the restricted
#'   maximum likelihood (REML) estimation and the other is
#'   \code{\link[climr]{fit_bayes}} by using the Bayesian inference with the
#'   Markov chain Monte Carlo (MCMC) techniques. We apply CMIP5 data into both
#'   of the inferences but for CMIP6, we just use the Bayesian inference.
#'
#'   For interpreting the results from both of functions(i.e.
#'   \code{\link[climr]{fit_reml}} and \code{\link[climr]{fit_bayes}}), there
#'   are commensurate S3 methods for \code{print}, \code{summary} and
#'   \code{print.summary}. We are setting result of \code{\link[climr]{fit_reml}}
#'   in class of \code{climr_reml} and result of \code{\link[climr]{fit_bayes}}
#'   in class of \code{climr_bayes}.
#
#' @references Northrop, P.J. and R.E. Chandler (2014). Quantifying
#'   Sources of Uncertainty in Projections of Future Climate.
#'   \emph{Journal of Climate}, \strong{27}, 8793-8808.
#'   \url{https://doi.org/10.1175/JCLI-D-14-00265.1}.
#'
#' @seealso \code{\link[climr]{cmip5_temp1}} and \code{\link[climr]{cmip5_temp2}}
#'  for CMIP5 data.
#' @seealso \code{\link[climr]{cmip6_pr}} and \code{\link[climr]{cmip6_tas}} for
#'   CMIP6 data.
#' @seealso \code{\link[climr]{plot.cmip5}} is a S3 plot method for CMIP5 data
#' whose class is \code{cmip5}.
#' @seealso \code{\link[climr]{plot.cmip6}} is a S3 plot method for CMIP6 data
#' whose class is \code{cmip6}.
#' @seealso \code{\link[climr]{fit_reml}} for the REML inference to obtain the
#' estimates.
#' @seealso \code{\link[climr]{print.climr_reml}} for \code{print} of result in
#' class of \code{climr_reml}.
#' @seealso \code{\link[climr]{summary.climr_reml}} for \code{summary} of result
#' in class of \code{climr_reml}.
#' @seealso \code{\link[climr]{print.summary.climr_reml}} for \code{print.summary}
#' of result in class of \code{climr_reml}.
#' @seealso \code{\link[climr]{fit_bayes}} for the Bayesian inference to obtain
#' the estimates.
#' @seealso \code{\link[climr]{print.climr_bayes}} for \code{print} of result
#' in class of \code{climr_bayes}.
#' @seealso \code{\link[climr]{summary.climr_bayes}} for \code{summary} of result
#' in class of \code{climr_bayes}.
#' @seealso \code{\link[climr]{print.summary.climr_bayes}} for \code{print.summary}
#' of result in class of \code{climr_bayes}.
#'
#' @docType package
#' @name climr
#' @importFrom graphics plot
NULL

#' Mid 21st Century Global Temperature Projection Data
#'
#' Indices of global temperature change from late 20th century (1970-1999)
#' to mid 21st century (2020-2049) based on data produced by the Fifth
#' Coupled Model Intercomparison Project (CMIP5).
#'
#' The data frame \code{cmip5_temp1} data frame has 270 rows and 4 columns.
#' Each row relates to a climate projection run from one of 38 different
#' General Circulation Models (GCMs) under a particular
#' Representative Concentration Pathway (RCP).
#' Use \code{table(cmip5_temp1[, c("GCM", "RCP")])} to see the numbers of
#' runs under each RCP for each GCM.
#' See Van Vuuren et al (2011) for an overview of RCPs
#' and Northrop and Chandler (2014) for analyses of a similar
#' older dataset (CMIP3).
#' Column 1 contains the anomaly of the mean global temperature over
#' the time period 2020-2049 relative to the mean global temperature
#' over 1970-1999, i.e. the latter subtracted from the former.
#' Column 2 contains an abbreviation for the name of the climate modelling
#' research group and the GCM.
#' Column 3 contains the RCP in the format \code{rcpxx} where \code{xx}
#' is a radiative forcing level resulting from an anticipated future
#' greenhouse gas emissions.
#' Column 4 is the simulation run number.
#' @format A data frame with 270 rows and 4 columns.
#'   \itemize{
#'     \item{Column 1, \code{index}: }{anomaly of 2020-2049 mean relative to
#'       the 1970-1999 mean.}
#'     \item{Column 2, \code{GCM}: }{Abbreviated name of General Circulation
#'       Model.}
#'     \item{Column 3, \code{RCP}: }{Representative Concentration Pathway.
#'       One of rcp26, rcp45, rcp60, rcp85.}
#'     \item{Column 4, \code{run}: }{Simulation run number.}
#'  }
#' @source The raw data from which the indices are calculated are monthly
#'   CMIP5 scenario runs for global surface air temperature (tas)
#'   downloaded from the KNMI Climate Explorer (\url{https://climexp.knmi.nl/})
#'   on 4/3/2015.
#' @references Northrop, P.J. and R.E. Chandler (2014). Quantifying
#'   Sources of Uncertainty in Projections of Future Climate.
#'   \emph{Journal of Climate}, \strong{27}, 8793-8808.
#'   \url{https://doi.org/10.1175/JCLI-D-14-00265.1}
#' @references Van Vuuren, D. P., Edmonds, J., Kainuma, M., Riahi, K.
#'   Thomson, A., Hibbard, K., Hurtt, G. C., Kram, T., Krey, V.,
#'   Lamarque, J.-F. (2011). The representative concentration pathways:
#'   an overview. \emph{Climatic change}, \strong{109}, 5-31.
#'   \url{https://doi.org/10.1007/s10584-011-0148-z}
"cmip5_temp1"

#' Late 21st Century Global Temperature Projection Data
#'
#' Indices of global temperature change from late 20th century (1970-1999)
#' to late 21st century (2069-2098) based on data produced by the Fifth
#' Coupled Model Intercomparison Project (CMIP5).
#'
#' The data frame \code{cmip5_temp2} data frame has 270 rows and 4 columns.
#' Each row relates to a climate projection run from one of 38 different
#' General Circulation Models (GCMs) under a particular
#' Representative Concentration Pathway (RCP).
#' Use \code{table(cmip5_temp2[, c("GCM", "RCP")])} to see the numbers of
#' runs under each RCP for each GCM.
#' See Van Vuuren et al (2011) for an overview of RCPs
#' and Northrop and Chandler (2014) for analyses of a similar
#' older dataset (CMIP3).
#' Column 1 contains the anomaly of the mean global temperature over
#' the time period 2069-2098 relative to the mean global temperature
#' over 1970-1999, i.e. the latter subtracted from the former.
#' Column 2 contains an abbreviation for the name of the climate modelling
#' research group and the GCM.
#' Column 3 contains the RCP in the format \code{rcpxx} where \code{xx}
#' is a radiative forcing level resulting from an anticipated future
#' greenhouse gas emissions.
#' Column 4 is the simulation run number.
#' @format A data frame with 270 rows and 4 columns.
#'   \itemize{
#'     \item{Column 1, \code{index}: }{anomaly of 2069-2098 mean relative to
#'       the 1970-1999 mean.}
#'     \item{Column 2, \code{GCM}: }{Abbreviated name of General Circulation
#'       Model.}
#'     \item{Column 3, \code{RCP}: }{Representative Concentration Pathway.
#'       One of rcp26, rcp45, rcp60, rcp85.}
#'     \item{Column 4, \code{run}: }{Simulation run number.}
#'  }
#' @source The raw data from which the indices are calculated are monthly
#'   CMIP5 scenario runs for global surface air temperature (tas)
#'   downloaded from the KNMI Climate Explorer (\url{https://climexp.knmi.nl/})
#'   on 4/3/2015.
#' @references Northrop, P.J. and R.E. Chandler (2014). Quantifying
#'   Sources of Uncertainty in Projections of Future Climate.
#'   \emph{Journal of Climate}, \strong{27}, 8793-8808.
#'   \url{https://doi.org/10.1175/JCLI-D-14-00265.1}
#' @references Van Vuuren, D. P., Edmonds, J., Kainuma, M., Riahi, K.
#'   Thomson, A., Hibbard, K., Hurtt, G. C., Kram, T., Krey, V.,
#'   Lamarque, J.-F. (2011). The representative concentration pathways:
#'   an overview. \emph{Climatic change}, \strong{109}, 5-31.
#'   \url{https://doi.org/10.1007/s10584-011-0148-z}
"cmip5_temp2"

#' CMIP6 21st Century Temperature Projection Data
#'
#' Indices of temperature change from late 20th century (1970-1999)
#' to mid 21st century (2020-2049) and to late 21st century (2070-2099) based
#' on data produced by the Sixth Coupled Model Intercomparison Project (CMIP6).
#' The climate variable is surface air temperature (tas).
#'
#' The data frame \code{cmip6_tas} data frame has 7268 rows and 6 columns.
#' Each row relates to a climate projection run from one of 11 different
#' General Circulation Models (GCMs) under a particular Shared Socio-Economic
#' Pathway (SSP) for a particular region of the world (and one of two time
#' horizons.
#'
#' See Northrop and Chandler (2014) for analyses of a similar older dataset
#' (CMIP3).
#' @format A data frame with 7268 rows and 6 columns.
#'   \itemize{
#'     \item{Column 1, \code{index}: }{anomaly of the 2020-2049 mean
#'       (if \code{period = 1} in column 6) or the 2070-2099 mean
#'       (if \code{period = 2} in column 6) relative to the 1970-1999 mean.
#'       That is, the latter subtracted from the former.}
#'     \item{Column 2, \code{GCM}: }{Abbreviated name of General Circulation
#'       Model.  There are 11 different GCMs.}
#'     \item{Column 3, \code{SSP}: }{Shared Socio-Economic Pathway.
#'       One of ssp119, ssp126, ssp245, ssp370, ssp434, ssp460, ssp534-over,
#'       ssp585.}
#'     \item{Column 4, \code{run}: }{Indicator of the simulation run number.}
#'     \item{Column 5, \code{region}: }{Region of the world: 0 for global,
#'       1-22 for regions.}
#'     \item{Column 6, \code{period}: }{Time horizon: 1 for mid 21st century
#'       (2020-2049); 2 for late 21st century (2070-2099).}
#'  }
#' @source The raw data from which the indices are calculated were downloaded
#'   via the World Climate Research Programme (WCRP) CMIP6 homepage
#'   \url{https://esgf-node.llnl.gov/projects/cmip6/}.
#' @references Northrop, P.J. and R.E. Chandler (2014). Quantifying
#'   Sources of Uncertainty in Projections of Future Climate.
#'   \emph{Journal of Climate}, \strong{27}, 8793-8808.
#'   \url{https://doi.org/10.1175/JCLI-D-14-00265.1}
"cmip6_tas"

#' CMIP6 21st Century Precipitation Projection Data
#'
#' Indices of percentage precipitation change from late 20th century (1970-1999)
#' to mid 21st century (2020-2049) and to late 21st century (2070-2099) based
#' on data produced by the Sixth Coupled Model Intercomparison Project (CMIP6).
#' The climate variable is precipitation flux (pr).
#'
#' The data frame \code{cmip6_pr} data frame has 7268 rows and 6 columns.
#' Each row relates to a climate projection run from one of 11 different
#' General Circulation Models (GCMs) under a particular Shared Socio-Economic
#' Pathway (SSP) for a particular region of the world (and one of two time
#' horizons.
#'
#' See Northrop and Chandler (2014) for analyses of a similar older dataset
#' (CMIP3).
#' @format A data frame with 7268 rows and 6 columns.
#'   \itemize{
#'     \item{Column 1, \code{index}: }{percentage anomaly of the 2020-2049 mean
#'       (if \code{period = 1} in column 6) or the 2070-2099 mean
#'       (if \code{period = 2} in column 6) relative to the 1970-1999 mean.
#'       That is, the latter subtracted from the former divided by former,
#'       expressed as a percentage.}
#'     \item{Column 2, \code{GCM}: }{Abbreviated name of General Circulation
#'       Model.  There are 11 different GCMs.}
#'     \item{Column 3, \code{SSP}: }{Shared Socio-Economic Pathway.
#'       One of ssp119, ssp126, ssp245, ssp370, ssp434, ssp460, ssp534-over,
#'       ssp585.}
#'     \item{Column 4, \code{run}: }{Indicator of the simulation run number.}
#'     \item{Column 5, \code{region}: }{Region of the world: 0 for global,
#'       1-22 for regions.}
#'     \item{Column 6, \code{period}: }{Time horizon: 1 for mid 21st century
#'       (2020-2049); 2 for late 21st century (2070-2099).}
#'  }
#' @source The raw data from which the indices are calculated were downloaded
#'   via the World Climate Research Programme (WCRP) CMIP6 homepage
#'   \url{https://esgf-node.llnl.gov/projects/cmip6/}.
#' @references Northrop, P.J. and R.E. Chandler (2014). Quantifying
#'   Sources of Uncertainty in Projections of Future Climate.
#'   \emph{Journal of Climate}, \strong{27}, 8793-8808.
#'   \url{https://doi.org/10.1175/JCLI-D-14-00265.1}
"cmip6_pr"
