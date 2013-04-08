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

  estimateDispersionsMAP(object, minDisp = 1e-08,
    quantilesForPrior, kappa_0 = 1, dispTol = 1e-04,
    maxit = 100)
}
\arguments{
  \item{object}{a DESeqSummarizedExperiment}

  \item{fitType}{either "parametric", "local", or "geoMean"
  for the type of fitting of dispersions to the mean
  intensity. See \code{\link{estimateDispersions}} for
  description.}

  \item{minDisp}{small value for the minimum dispersion, to
  allow for calculations in log scale, one decade above
  this value is used as a test for inclusion in
  mean-dispersion fitting}

  \item{quantilesForPrior}{numeric vector, e.g.
  \code{c(lower, upper)} specifying by quantiles of row
  mean which rows to use for calculating the prior width.
  If not provided, for fitType = "parametric" or "local",
  the quantiles \code{c(.5,1)} are used, otherwise
  \code{c(0,1)}.}

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
  a DESeqSummarizedExperiment with gene-wise, fitted, or
  final MAP dispersion estimates in the metadata columns of
  the object.
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
dse <- makeExampleDESeqSummarizedExperiment()
dse <- estimateSizeFactors(dse)
dse <- estimateDispersionsGeneEst(dse)
mcols(dse)$dispFit <- rep(median(mcols(dse)$dispGeneEst,na.rm=TRUE),nrow(dse))
dse <- estimateDispersionsMAP(dse, quantilesForPrior=c(0,1))
plotDispEsts(dse)
}
\seealso{
  \code{\link{estimateDispersions}}
}

