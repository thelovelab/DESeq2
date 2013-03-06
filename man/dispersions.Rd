\docType{methods}
\name{dispersions}
\alias{dispersions}
\alias{dispersions,DESeqSummarizedExperiment-method}
\alias{dispersions<-}
\alias{dispersions<-,DESeqSummarizedExperiment,numeric-method}
\title{Accessor functions for the dispersion estimates in a DESeqSummarizedExperiment
object.}
\usage{
  \S4method{dispersions}{DESeqSummarizedExperiment}(object)

  \S4method{dispersions}{DESeqSummarizedExperiment,numeric}(object)<-value
}
\arguments{
  \item{object}{a \code{DESeqSummarizedExperiment} object.}

  \item{value}{the dispersions to use for the negative
  binomial modeling}
}
\description{
  The dispersions for each row of the
  DESeqSummarizedExperiment.  Generally, these should be
  set only by \code{\link{estimateDispersions}}.
}
\examples{
example("estimateDispersions")
dispersions(dse)
}
\author{
  Simon Anders
}
\seealso{
  \code{\link{estimateDispersions}}
}

