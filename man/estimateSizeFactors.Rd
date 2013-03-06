\docType{methods}
\name{estimateSizeFactors}
\alias{estimateSizeFactors}
\alias{estimateSizeFactors,DESeqSummarizedExperiment-method}
\title{Estimate the size factors for a DESeqSummarizedExperiment}
\usage{
  \S4method{estimateSizeFactors}{DESeqSummarizedExperiment}(object,locfunc=median)
}
\arguments{
  \item{object}{a DESeqSummarizedExperiment}

  \item{locfunc}{a function to compute a location for a
  sample. By default, the median is used. However,
  especially for low counts, the
  \code{\link[genefilter]{shorth}} function from the
  genefilter package may give better results.}
}
\value{
  The DESeqSummarizedExperiment passed as parameters, with
  the size factors filled in.
}
\description{
  Estimate the size factors for a DESeqSummarizedExperiment
}
\details{
  This function estimates the size factors and stores the
  information which can be accessed using
  \code{\link{sizeFactors}}

  Typically, the function is called with the idiom:

  \code{dse <- estimateSizeFactors(dse)}

  See \code{\link{DESeq}} for a description of the use of
  size factors in the GLM. You need to call this function
  after \code{\link{DESeqSummarizedExperiment}} unless you
  have manually specified \code{\link{sizeFactors}}.
  Alternatively, gene-specific normalization factors for
  each sample can be provided using
  \code{\link{normalizationFactors}} which will always
  preempt \code{\link{sizeFactors}} in calculations.

  Internally, the function calls
  \code{\link{estimateSizeFactorsForMatrix}}, which
  provides more details on the calculation.
}
\examples{
dse <- makeExampleDESeqSummarizedExperiment()
dse <- estimateSizeFactors( dse )
sizeFactors( dse )
}
\author{
  Simon Anders
}
\seealso{
  \code{\link{estimateSizeFactorsForMatrix}}
}

