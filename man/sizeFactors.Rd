\docType{methods}
\name{sizeFactors}
\alias{sizeFactors}
\alias{sizeFactors,DESeqDataSet-method}
\alias{sizeFactors<-,DESeqDataSet,numeric-method}
\title{Accessor functions for the 'sizeFactors' information in a DESeqDataSet
object.}
\usage{
  \S4method{sizeFactors}{DESeqDataSet}(object)

  \S4method{sizeFactors}{DESeqDataSet,numeric}(object)<-value
}
\arguments{
  \item{object}{a \code{DESeqDataSet} object.}

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
dds <- makeExampleDESeqDataSet()
dds <- estimateSizeFactors( dds )
sizeFactors(dds)
}
\author{
  Simon Anders
}
\seealso{
  \code{\link{estimateSizeFactors}}
}

