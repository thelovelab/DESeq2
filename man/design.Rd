\docType{methods}
\name{design}
\alias{design}
\alias{design,DESeqSummarizedExperiment-method}
\alias{design<-,DESeqSummarizedExperiment,formula-method}
\title{Accessors for the 'design' slot of a DESeqSummarizedExperiment object.}
\usage{
  \S4method{design}{DESeqSummarizedExperiment}(object)

  \S4method{design}{DESeqSummarizedExperiment,formula}(object)<-value
}
\arguments{
  \item{object}{a \code{DESeqSummarizedExperiment} object}

  \item{value}{a \code{formula} used for estimating
  dispersion and fitting negative binomial GLMs}
}
\description{
  Accessors for the 'design' slot of a
  DESeqSummarizedExperiment object.
}
\examples{
dse <- makeExampleDESeqSummarizedExperiment()
design(dse) <- formula(~ 1)
}

