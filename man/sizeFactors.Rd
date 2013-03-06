\docType{methods}
\name{sizeFactors}
\alias{sizeFactors}
\alias{sizeFactors,DESeqSummarizedExperiment-method}
\alias{sizeFactors<-,DESeqSummarizedExperiment,numeric-method}
\title{Accessor functions for the 'sizeFactors' information in a DESeqSummarizedExperiment
object.}
\usage{
  \S4method{sizeFactors}{DESeqSummarizedExperiment}(object)

  \S4method{sizeFactors}{DESeqSummarizedExperiment,numeric}(object)<-value
}
\arguments{
  \item{object}{a \code{DESeqSummarizedExperiment} object.}

  \item{value}{a numeric vector, one size factor for each
  column in the count data.}
}
\description{
  The sizeFactors vector assigns to each column of the
  count matrix a value, the size factor, such that count
  values in the columns can be brought to a common scale by
  dividing by the corresponding size factor.  See
  \code{\link{DESeq}} for a description of the use of size
  factors. If gene-specific normalization is desired for
  each sample, use \code{\link{normalizationFactors}}.
}
\examples{
dse <- makeExampleDESeqSummarizedExperiment()
dse <- estimateSizeFactors( dse )
sizeFactors(dse)
}
\author{
  Simon Anders
}
\seealso{
  \code{\link{estimateSizeFactors}}
}

