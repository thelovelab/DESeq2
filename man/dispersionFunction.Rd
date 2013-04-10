\docType{methods}
\name{dispersionFunction}
\alias{dispersionFunction}
\alias{dispersionFunction,DESeqDataSet-method}
\alias{dispersionFunction<-}
\alias{dispersionFunction<-,DESeqDataSet,function-method}
\title{Accessors for the 'dispersionFunction' slot of a DESeqDataSet object.}
\usage{
  \S4method{dispersionFunction}{DESeqDataSet}(object)

  \S4method{dispersionFunction}{DESeqDataSet,function}(object)<-value
}
\arguments{
  \item{object}{a \code{DESeqDataSet} object.}

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
dispersionFunction(dds)
}

