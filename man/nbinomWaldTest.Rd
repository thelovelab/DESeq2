\name{nbinomWaldTest}
\alias{nbinomWaldTest}
\title{Wald test for the GLM coefficients}
\usage{
nbinomWaldTest(object, betaPrior = TRUE, betaPriorVar, modelMatrixType,
  maxit = 100, useOptim = TRUE, quiet = FALSE, useT = FALSE, df,
  useQR = TRUE)
}
\arguments{
  \item{object}{a DESeqDataSet}

  \item{betaPrior}{whether or not to put a zero-mean normal
  prior on the non-intercept coefficients (Tikhonov/ridge
  regularization)}

  \item{betaPriorVar}{a vector with length equal to the
  number of model terms including the intercept. which if
  missing is estimated from the rows which do not have any
  zeros}

  \item{modelMatrixType}{either "standard" or "expanded",
  which describe how the model matrix, X of the formula in
  \code{\link{DESeq}}, is formed. "standard" is as created
  by \code{model.matrix} using the design formula.
  "expanded" includes an indicator variable for each level
  of factors with 3 or more levels, in order to ensure
  symmetric behavior of the prior on log2 fold changes.
  betaPrior must be set to TRUE in order for expanded model
  matrices to be fit.}

  \item{maxit}{the maximum number of iterations to allow
  for convergence of the coefficient vector}

  \item{useOptim}{whether to use the native optim function
  on rows which do not converge within maxit}

  \item{quiet}{whether to print messages at each step}

  \item{useT}{whether to use a t-distribution as a null
  distribution, for significance testing of the Wald
  statistics. If FALSE, a standard normal null distribution
  is used.}

  \item{df}{the degrees of freedom for the t-distribution}

  \item{useQR}{whether to use the QR decomposition on the
  design matrix X while fitting the GLM}
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
likelihood estimates for GLM coefficients are calculated; a
zero-mean normal prior distribution is assumed; the
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
normal distribution. From examination of Wald statistics
for real datasets, the effect of the prior on dispersion
estimates results in a Wald statistic distribution which is
approximately normal.

When interaction terms are present, the prior on log fold
changes the calculated beta prior variance will only be
used for the interaction terms (non-interaction terms
receive a wide prior variance of 1e6). In the case of
interaction terms and factors with 3 or more levels present
in the design formula, a moderately wide prior variance of
1e3 will be used on non-interaction terms, to allow for
convergence to the maximum of the posterior.

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
\code{\link{DESeq}}, \code{\link{nbinomLRT}}
}

