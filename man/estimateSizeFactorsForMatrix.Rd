\name{estimateSizeFactorsForMatrix}
\alias{estimateSizeFactorsForMatrix}
\title{Low-level function to estimate size factors with robust regression.}
\usage{
  estimateSizeFactorsForMatrix(counts, locfunc = median)
}
\arguments{
  \item{counts}{a matrix or data frame of counts, i.e.,
  non-negative integer values}

  \item{locfunc}{a function to compute a location for a
  sample. By default, the median is used. However,
  especially for low counts, the
  \code{\link[genefilter]{shorth}} function from genefilter
  may give better results.}
}
\value{
  a vector with the estimates size factors, one element per
  column
}
\description{
  Given a matrix or data frame of count data, this function
  estimates the size factors as follows: Each column is
  divided by the geometric means of the rows. The median
  (or, if requested, another location estimator) of these
  ratios (skipping the genes with a geometric mean of zero)
  is used as the size factor for this column. Typically,
  you will not call this function directly, but use
  \code{\link{estimateSizeFactors}}.
}
\examples{
dds <- makeExampleDESeqDataSet()
estimateSizeFactorsForMatrix(counts(dds))
}
\author{
  Simon Anders
}
\seealso{
  \code{\link{estimateSizeFactors}}
}

