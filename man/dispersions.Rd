\docType{methods}
\name{dispersions}
\alias{dispersions}
\alias{dispersions,DESeqDataSet-method}
\alias{dispersions<-}
\alias{dispersions<-,DESeqDataSet,numeric-method}
\title{Accessor functions for the dispersion estimates in a DESeqDataSet
object.}
\usage{
  \S4method{dispersions}{DESeqDataSet}(object)

  \S4method{dispersions}{DESeqDataSet,numeric}(object)<-value
}
\arguments{
  \item{object}{a \code{DESeqDataSet} object.}

  \item{value}{the dispersions to use for the negative
  binomial modeling}
}
\description{
  The dispersions for each row of the DESeqDataSet.
  Generally, these should be set only by
  \code{\link{estimateDispersions}}.
}
\examples{
example("estimateDispersions")
dispersions(dds)
}
\author{
  Simon Anders
}
\seealso{
  \code{\link{estimateDispersions}}
}

