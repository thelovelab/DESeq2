\name{nbinomLRT}
\alias{nbinomLRT}
\title{Likelihood ratio test (chi-squared test) for GLMs}
\usage{
  nbinomLRT(object, full = design(object), reduced,
    pAdjustMethod = "BH", maxit = 100, useOptim = TRUE,
    quiet = FALSE)
}
\arguments{
  \item{object}{a DESeqDataSet}

  \item{full}{the full model formula, this should be the
  formula in \code{design(object)}}

  \item{reduced}{a reduced formula to compare against, e.g.
  the full model with a variable of interest removed}

  \item{pAdjustMethod}{the method to use for adjusting
  p-values, see \code{?p.adjust}}

  \item{maxit}{the maximum number of iterations to allow
  for convergence of the coefficient vector}

  \item{useOptim}{whether to use the native optim function
  on rows which do not converge within maxit}

  \item{quiet}{whether to print messages at each step}
}
\value{
  a DESeqDataSet with new results columns accessible with
  the \code{\link{results}} function.  The coefficients and
  standard errors are reported on a log2 scale.
}
\description{
  This function tests for significance of change in
  deviance between a full and reduced model which are
  provided as \code{formula}. Fitting uses previously
  calculated \code{\link{sizeFactors}} (or
  \code{\link{normalizationFactors}}) and dispersion
  estimates.
}
\details{
  The difference in deviance is compared to a chi-squared
  distribution with df = (reduced residual degrees of
  freedom - full residual degrees of freedom). This
  function is comparable to the \code{nbinomGLMTest} of the
  previous version of DESeq and an alternative to the
  default \code{\link{nbinomWaldTest}}.
}
\examples{
dds <- makeExampleDESeqDataSet()
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomLRT(dds, reduced = ~ 1)
res <- results(dds)
}
\seealso{
  \code{\link{nbinomWaldTest}}
}

