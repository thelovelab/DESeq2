\docType{methods}
\name{design}
\alias{design}
\alias{design,DESeqDataSet-method}
\alias{design<-,DESeqDataSet,formula-method}
\title{Accessors for the 'design' slot of a DESeqDataSet object.}
\usage{
  \S4method{design}{DESeqDataSet}(object)

  \S4method{design}{DESeqDataSet,formula}(object)<-value
}
\arguments{
  \item{object}{a \code{DESeqDataSet} object}

  \item{value}{a \code{formula} used for estimating
  dispersion and fitting negative binomial GLMs}
}
\description{
  Accessors for the 'design' slot of a DESeqDataSet object.
}
\examples{
dds <- makeExampleDESeqDataSet()
design(dds) <- formula(~ 1)
}

