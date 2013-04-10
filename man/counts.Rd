\docType{methods}
\name{counts}
\alias{counts}
\alias{counts,DESeqDataSet-method}
\alias{counts<-,DESeqDataSet,matrix-method}
\title{Accessors for the 'counts' slot of a DESeqDataSet object.}
\usage{
  \S4method{counts}{DESeqDataSet}(object,normalized=FALSE)

  \S4method{counts}{DESeqDataSet,matrix}(object)<-value
}
\arguments{
  \item{object}{a \code{DESeqDataSet} object.}

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
dds <- makeExampleDESeqDataSet()
head(counts(dds))
}
\author{
  Simon Anders
}
\seealso{
  \code{\link{sizeFactors}},
  \code{\link{normalizationFactors}}
}

