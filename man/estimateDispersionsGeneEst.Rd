\name{estimateDispersionsGeneEst}
\alias{estimateDispersionsFit}
\alias{estimateDispersionsGeneEst}
\alias{estimateDispersionsMAP}
\title{Low-level functions to fit dispersion estimates}
\usage{
  estimateDispersionsGeneEst(object, minDisp = 1e-08,
    kappa_0 = 1, dispTol = 1e-06, maxit = 100)

  estimateDispersionsFit(object,
    fitType = c("parametric", "local", "geoMean"),
    minDisp = 1e-08)

  estimateDispersionsMAP(object, outlierSD = 2, priorVar,
    minDisp = 1e-08, kappa_0 = 1, dispTol = 1e-06,
    maxit = 100)
}
\arguments{
  \item{object}{a DESeqDataSet}

  \item{fitType}{either "parametric", "local", or "geoMean"
  for the type of fitting of dispersions to the mean
  intensity. See \code{\link{estimateDispersions}} for
  description.}

  \item{outlierSD}{the number of prior standard deviations
  above the prior mean (fitted value), above which
  dispersion estimates will be labelled outliers. Outliers
  will keep their original value and not be shrunk using
  the prior.}

  \item{priorVar}{the variance of the normal prior on the
  log dispersions. If not supplied, this is calculated as
  the difference between the mean squared residuals of
  gene-wise estimates to the fitted dispersion and the
  expected sampling variance of the log dispersion}

  \item{minDisp}{small value for the minimum dispersion, to
  allow for calculations in log scale, one decade above
  this value is used as a test for inclusion in
  mean-dispersion fitting}

  \item{kappa_0}{control parameter used in setting the
  initial proposal in backtracking search, higher kappa_0
  results in larger steps}

  \item{dispTol}{control parameter to test for convergence
  of log dispersion, stop when increase in log posterior is
  less than dispTol}

  \item{maxit}{control parameter: maximum number of
  iterations to allow for convergence}
}
\value{
  a DESeqDataSet with gene-wise, fitted, or final MAP
  dispersion estimates in the metadata columns of the
  object.
}
\description{
  Normal users should instead use
  \code{\link{estimateDispersions}}. These low-level
  functions are called by
  \code{\link{estimateDispersions}}, but are exported and
  documented for non-standard usage. For instance, it is
  possible to replace fitted values with a custom fit and
  continue with the maximum a posteriori dispersion
  estimation, as demonstrated in the examples below.
}
\examples{
dds <- makeExampleDESeqDataSet()
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
mcols(dds)$dispFit <- rep(median(mcols(dds)$dispGeneEst,na.rm=TRUE),nrow(dds))
dds <- estimateDispersionsMAP(dds)
plotDispEsts(dds)
}
\seealso{
  \code{\link{estimateDispersions}}
}

