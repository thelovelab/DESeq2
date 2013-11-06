\name{replaceOutliersWithTrimmedMean}
\alias{replaceOutliersWithTrimmedMean}
\title{Replace outliers with trimmed mean}
\usage{
  replaceOutliersWithTrimmedMean(dds, trim = 0.2,
    cooksCutoff)
}
\arguments{
  \item{dds}{a DESeqDataSet object, which has already been
  processed by either DESeq, nbinomWaldTest or nbinomLRT,
  and therefore contains a matrix 'cooks' contained in
  assays(dds). These are the Cook's distances which will be
  used to define outlier counts.}

  \item{trim}{the fraction (0 to 0.5) of observations to be
  trimmed from each end of the normalized counts for a gene
  before the mean is computed}

  \item{cooksCutoff}{the threshold for defining an outlier
  to be replaced. Defaults to the .99 quantile of the F(p,
  m - p) distribution, where p is the number of parameters
  and m is the number of samples.}
}
\description{
  The \code{\link{DESeq}} function calculates a diagnostic
  measure called Cook's distance for every gene and every
  sample. The \code{\link{results}} function then sets the
  p-values to \code{NA} for genes which contain an outlying
  count as defined by a Cook's distance above a threshold.
  With many degrees of freedom, i.e. many more samples than
  number of parameters to be estimated-- it might be
  undesirable to remove entire genes from the analysis just
  because their data include a single count outlier. An
  alternate strategy is to replace the outlier counts with
  the trimmed mean over all samples, adjusted by the size
  factor or normalization factor for that sample. The
  following simple function performs this replacement for
  the user. For more information on Cook's distance, please
  see the two sections of the vignette: 'Dealing with count
  outliers' and 'Count outlier detection'.
}
\examples{
dds <- makeExampleDESeqDataSet(n=100)
dds <- DESeq(dds)
ddsReplace <- replaceOutliersWithTrimmedMean(dds)
}
\seealso{
  \code{\link{DESeq}}
}

