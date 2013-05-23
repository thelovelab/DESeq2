\name{nbinomWaldTest}
\alias{nbinomWaldTest}
\title{Wald test for the GLM coefficients}
\usage{
  nbinomWaldTest(object, betaPrior = TRUE,
    pAdjustMethod = "BH", priorSigmasq, maxit = 100,
    useOptim = TRUE, quiet = FALSE, useT = FALSE, dfCoef)
}
\arguments{
  \item{object}{a DESeqDataSet}

  \item{betaPrior}{whether or not to put a zero-mean normal
  prior on the non-intercept coefficients (Tikhonov/ridge
  regularization)}

  \item{pAdjustMethod}{the method to use for adjusting
  p-values, see \code{?p.adjust}}

  \item{priorSigmasq}{a vector with length equal to the
  number of model terms including the intercept. which if
  missing is estimated from the rows which do not have any
  zeros}

  \item{maxit}{the maximum number of iterations to allow
  for convergence of the coefficient vector}

  \item{useOptim}{whether to use the native optim function
  on rows which do not converge within maxit}

  \item{quiet}{whether to print messages at each step}

  \item{useT}{whether to use as a null distribution, for
  significance testing of the Wald statistics, a
  t-distribution with degrees of freedom calculated with
  the dispersion prior variance. If FALSE, a standard
  normal null distribution is used.}

  \item{dfCoef}{a parameter used, if useT is set to TRUE,
  to connect the dispersion prior variance to the degrees
  of freedom of the t-distribution for Wald statistics: DF
  = dfCoef/priorVar + (m - p)}
}
\value{
  a DESeqDataSet with results columns accessible with the
  \code{\link{results}} function.  The coefficients and
  standard errors are reported on a log2 scale.
}
\description{
  This function tests for significance of coefficients in a
  negative binomial GLM, using previously calculated
  \code{\link{sizeFactors}} (or
  \code{\link{normalizationFactors}}) and dispersion
  estimates.  See \code{\link{DESeq}} for the GLM formula.
}
\details{
  The fitting proceeds as follows: standard maximum
  likelihood estimates for GLM coefficients are calculated;
  a zero-mean normal prior distribution is assumed; the
  variance of the prior distribution for each non-intercept
  coefficient is calculated as the mean squared maximum
  likelihood estimates over the genes which do not contain
  zeros for some condition; the final coefficients are then
  maximum a posteriori estimates (using Tikhonov/ridge
  regularization) using this prior. The use of a prior has
  little effect on genes with high counts and helps to
  moderate the large spread in coefficients for genes with
  low counts.

  For calculating Wald test p-values, the coefficients are
  scaled by their standard errors and then compared to a
  normal distribution.  From examination of Wald statistics
  for real datasets, the effect of the prior on dispersion
  estimates results in a Wald statistic distribution which
  is approximately normal.

  The Wald test can be replaced with the
  \code{\link{nbinomLRT}} for an alternative test of
  significance.
}
\examples{
dds <- makeExampleDESeqDataSet()
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)
}
\seealso{
  \code{\link{nbinomLRT}}
}

