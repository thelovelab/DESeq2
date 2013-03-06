\docType{methods}
\name{dispersionFunction}
\alias{dispersionFunction}
\alias{dispersionFunction,DESeqSummarizedExperiment-method}
\alias{dispersionFunction<-}
\alias{dispersionFunction<-,DESeqSummarizedExperiment,function-method}
\title{Accessors for the 'dispersionFunction' slot of a DESeqSummarizedExperiment object.}
\usage{
  \S4method{dispersionFunction}{DESeqSummarizedExperiment}(object)

  \S4method{dispersionFunction}{DESeqSummarizedExperiment,function}(object)<-value
}
\arguments{
  \item{object}{a \code{DESeqSummarizedExperiment} object.}

  \item{value}{a \code{function}}
}
\description{
  The dispersion function is calculated by
  \code{\link{estimateDispersions}} and used by
  \code{\link{varianceStabilizingTransformation}}.
  Parametric dispersion fits store the coefficients of the
  fit as attributes in this slot.
}
\examples{
example("estimateDispersions")
dispersionFunction(dse)
}

