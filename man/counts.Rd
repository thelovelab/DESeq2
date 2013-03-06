\docType{methods}
\name{counts}
\alias{counts}
\alias{counts,DESeqSummarizedExperiment-method}
\alias{counts<-,DESeqSummarizedExperiment,matrix-method}
\title{Accessors for the 'counts' slot of a DESeqSummarizedExperiment object.}
\usage{
  \S4method{counts}{DESeqSummarizedExperiment}(object,normalized=FALSE)

  \S4method{counts}{DESeqSummarizedExperiment,matrix}(object)<-value
}
\arguments{
  \item{object}{a \code{DESeqSummarizedExperiment} object.}

  \item{normalized}{logical indicating whether or not to
  divide the counts by the size factors or normalization
  factors before returning (normalization factors always
  preempt size factors)}

  \item{value}{an integer matrix}
}
\description{
  The counts slot holds the count data as a matrix of
  non-negative integer count values, one row for each
  observational unit (gene or the like), and one column for
  each sample.
}
\examples{
dse <- makeExampleDESeqSummarizedExperiment()
head(counts(dse))
}
\author{
  Simon Anders
}
\seealso{
  \code{\link{sizeFactors}},
  \code{\link{normalizationFactors}}
}

