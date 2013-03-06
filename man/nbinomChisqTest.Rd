\name{nbinomChisqTest}
\alias{nbinomChisqTest}
\title{Chi-squared test for GLMs}
\usage{
  nbinomChisqTest(object, full = design(object), reduced,
    pAdjustMethod = "BH")
}
\arguments{
  \item{object}{a DESeqSummarizedExperiment}

  \item{full}{the full model formula, this should be the
  formula in \code{design(object)}}

  \item{reduced}{a reduced formula to compare against, e.g.
  the full model with a variable of interest removed}

  \item{pAdjustMethod}{the method to use for adjusting
  p-values, see \code{?p.adjust}}
}
\value{
  a DESeqSummarizedExperiment with new results columns
  accessible with the \code{\link{results}} function.  The
  coefficients and standard errors are reported on a log2
  scale.
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
dse <- makeExampleDESeqSummarizedExperiment()
dse <- estimateSizeFactors(dse)
dse <- estimateDispersions(dse)
dse <- nbinomChisqTest(dse, reduced = ~ 1)
res <- results(dse)
}
\seealso{
  \code{\link{nbinomWaldTest}}
}

