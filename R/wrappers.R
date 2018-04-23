# Fit dispersions for negative binomial GLM
#
# This function estimates the dispersion parameter (alpha) for negative binomial
# generalized linear models. The fitting is performed on the log scale.
#
# ySEXP n by m matrix of counts
# xSEXP m by k design matrix
# mu_hatSEXP n by m matrix, the expected mean values, given beta-hat
# log_alphaSEXP n length vector of initial guesses for log(alpha)
# log_alpha_prior_meanSEXP n length vector of the fitted values for log(alpha)
# log_alpha_prior_sigmasqSEXP a single numeric value for the variance of the prior
# min_log_alphaSEXP the minimum value of log alpha
# kappa_0SEXP a parameter used in calculting the initial proposal
#   for the backtracking search
#   initial proposal = log(alpha) + kappa_0 * deriv. of log lik. w.r.t. log(alpha)
# tolSEXP tolerance for convergence in estimates
# maxitSEXP maximum number of iterations
# usePriorSEXP boolean variable, whether to use a prior or just calculate the MLE
# weightsSEXP n by m matrix of weights
# useWeightsSEXP whether to use weights
#
# return a list with elements: log_alpha, iter, iter_accept, last_change, initial_lp, intial_dlp, last_lp, last_dlp, last_d2lp
fitDispWrapper <- function (ySEXP, xSEXP, mu_hatSEXP, log_alphaSEXP, log_alpha_prior_meanSEXP,
                            log_alpha_prior_sigmasqSEXP, min_log_alphaSEXP, kappa_0SEXP,
                            tolSEXP, maxitSEXP, usePriorSEXP, weightsSEXP, useWeightsSEXP) {
  # test for any NAs in arguments
  arg.names <- names(formals(fitDispWrapper))
  na.test <- sapply(mget(arg.names), function(x) any(is.na(x)))
  if (any(na.test)) stop(paste("in call to fitDisp, the following arguments contain NA:",
                               paste(arg.names[na.test],collapse=", ")))
  fitDisp(ySEXP=ySEXP, xSEXP=xSEXP, mu_hatSEXP=mu_hatSEXP,
          log_alphaSEXP=log_alphaSEXP, log_alpha_prior_meanSEXP=log_alpha_prior_meanSEXP,
          log_alpha_prior_sigmasqSEXP=log_alpha_prior_sigmasqSEXP,
          min_log_alphaSEXP=min_log_alphaSEXP, kappa_0SEXP=kappa_0SEXP,
          tolSEXP=tolSEXP, maxitSEXP=maxitSEXP, usePriorSEXP=usePriorSEXP,
          weightsSEXP=weightsSEXP, useWeightsSEXP=useWeightsSEXP)
}

# Fit dispersions by evaluating over grid
#
# This function estimates the dispersion parameter (alpha) for negative binomial
# generalized linear models. The fitting is performed on the log scale.
#
# ySEXP n by m matrix of counts
# xSEXP m by k design matrix
# mu_hatSEXP n by m matrix, the expected mean values, given beta-hat
# disp_gridSEXP the grid over which to estimate
# log_alpha_prior_meanSEXP n length vector of the fitted values for log(alpha)
# log_alpha_prior_sigmasqSEXP a single numeric value for the variance of the prior
# usePriorSEXP boolean variable, whether to use a prior or just calculate the MLE
# weightsSEXP n by m matrix of weights
# useWeightsSEXP whether to use weights
#
# return a list with elements: 
fitDispGridWrapper <- function(y, x, mu, logAlphaPriorMean, logAlphaPriorSigmaSq, usePrior,
                               weightsSEXP, useWeightsSEXP) {
  # test for any NAs in arguments
  arg.names <- names(formals(fitDispGridWrapper))
  na.test <- sapply(mget(arg.names), function(x) any(is.na(x)))
  if (any(na.test)) stop(paste("in call to fitDispGridWrapper, the following arguments contain NA:",
                               paste(arg.names[na.test],collapse=", ")))
  minLogAlpha <- log(1e-8)
  maxLogAlpha <- log(max(10, ncol(y)))
  dispGrid <- seq(from=minLogAlpha, to=maxLogAlpha, length=20)
  logAlpha <- fitDispGrid(ySEXP=y, xSEXP=x, mu_hatSEXP=mu, disp_gridSEXP=dispGrid,
                          log_alpha_prior_meanSEXP=logAlphaPriorMean,
                          log_alpha_prior_sigmasqSEXP=logAlphaPriorSigmaSq,
                          usePriorSEXP=usePrior,
                          weightsSEXP=weightsSEXP, useWeightsSEXP=useWeightsSEXP)$log_alpha
  exp(logAlpha)
}

# Fit beta coefficients for negative binomial GLM
#
# This function estimates the coefficients (betas) for negative binomial generalized linear models.
#
# ySEXP n by m matrix of counts
# xSEXP m by k design matrix
# nfSEXP n by m matrix of normalization factors
# alpha_hatSEXP n length vector of the disperion estimates
# contrastSEXP a k length vector for a possible contrast
# beta_matSEXP n by k matrix of the initial estimates for the betas
# lambdaSEXP k length vector of the ridge values
# weightsSEXP n by m matrix of weights
# useWeightsSEXP whether to use weights
# tolSEXP tolerance for convergence in estimates
# maxitSEXP maximum number of iterations
# useQRSEXP whether to use QR decomposition
#
# Note: at this level the betas are on the natural log scale
fitBetaWrapper <- function (ySEXP, xSEXP, nfSEXP, alpha_hatSEXP, contrastSEXP,
                            beta_matSEXP, lambdaSEXP, weightsSEXP, useWeightsSEXP,
                            tolSEXP, maxitSEXP, useQRSEXP, minmuSEXP) {
  if ( missing(contrastSEXP) ) {
    # contrast is not required, just give 1,0,0,...
    contrastSEXP <- c(1,rep(0,ncol(xSEXP)-1))
  }
  # test for any NAs in arguments
  arg.names <- names(formals(fitBetaWrapper))
  na.test <- sapply(mget(arg.names), function(x) any(is.na(x)))
  if (any(na.test)) stop(paste("in call to fitBeta, the following arguments contain NA:",
                               paste(arg.names[na.test],collapse=", ")))
  
  fitBeta(ySEXP=ySEXP, xSEXP=xSEXP, nfSEXP=nfSEXP, alpha_hatSEXP=alpha_hatSEXP,
          contrastSEXP=contrastSEXP, beta_matSEXP=beta_matSEXP,
          lambdaSEXP=lambdaSEXP, weightsSEXP=weightsSEXP, useWeightsSEXP=useWeightsSEXP,
          tolSEXP=tolSEXP, maxitSEXP=maxitSEXP, useQRSEXP=useQRSEXP, minmuSEXP=minmuSEXP)
}
