\name{varianceStabilizingTransformation}
\alias{getVarianceStabilizedData}
\alias{varianceStabilizingTransformation}
\title{Apply a variance stabilizing transformation (VST) to the count data}
\usage{
  varianceStabilizingTransformation(object, blind = TRUE)

  getVarianceStabilizedData(object)
}
\arguments{
  \item{object}{a DESeqDataSet, with \code{design(object)
  <- formula(~ 1)} and size factors (or normalization
  factors) and dispersions estimated using local or
  parametric \code{fitType}.}

  \item{blind}{logical, whether to blind the transformation
  to the experimental design. blind=TRUE should be used for
  comparing samples in an manner unbiased by prior
  information on samples, for example to perform sample QA
  (quality assurance). blind=FALSE should be used for
  transforming data for downstream analysis, where the full
  use of the design information should be made.}
}
\value{
  for \code{varianceStabilizingTransformation}, a
  \code{SummarizedExperiment}. for
  \code{getVarianceStabilizedData}, a \code{matrix} of the
  same dimension as the count data, containing the
  transformed values.
}
\description{
  This function calculates a variance stabilizing
  transformation (VST) from the fitted dispersion-mean
  relation(s) and then transforms the count data
  (normalized by division by the size factors or
  normalization factors), yielding a matrix of values which
  are now approximately homoskedastic (having constant
  variance along the range of mean values). The
  \code{\link{rlogTransformation}} is less sensitive to
  size factors, which can be an issue when size factors
  vary widely. This transformation is useful when checking
  for outliers or as input for machine learning techniques
  such as clustering or linear discriminant analysis.
}
\details{
  For each sample (i.e., column of \code{counts(dds)}), the
  full variance function is calculated from the raw
  variance (by scaling according to the size factor and
  adding the shot noise). We recommend an unsupervised
  estimation of the variance function, i.e., one ignoring
  conditions. This is performed by default, and can be
  modified using the 'unsupervised' argument.

  A typical workflow is shown in Section \emph{Variance
  stabilizing transformation} in the package vignette.

  If \code{\link{estimateDispersions}} was called with
  \code{fitType="parametric"}, a closed-form expression for
  the variance stabilizing transformation is used on the
  normalized count data. The expression can be found in the
  file \file{vst.pdf} which is distributed with the
  vignette.

  If \code{\link{estimateDispersions}} was called with
  \code{fitType="local"}, the reciprocal of the square root
  of the variance of the normalized counts, as derived from
  the dispersion fit, is then numerically integrated, and
  the integral (approximated by a spline function) is
  evaluated for each count value in the column, yielding a
  transformed value.

  In both cases, the transformation is scaled such that for
  large counts, it becomes asymptotically (for large
  values) equal to the logarithm to base 2.

  Limitations: In order to preserve normalization, the same
  transformation has to be used for all samples. This
  results in the variance stabilizition to be only
  approximate. The more the size factors differ, the more
  residual dependence of the variance on the mean you will
  find in the transformed data. As shown in the vignette,
  you can use the function \code{meanSdPlot} from the
  package \pkg{vsn} to see whether this is a problem for
  your data.
}
\examples{
dds <- makeExampleDESeqDataSet()
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
par(mfrow=c(1,2))
plot(rank(rowMeans(counts(dds))), genefilter::rowVars(log2(counts(dds)+1)), main="log2(x+1) transform")
plot(rank(rowMeans(assay(vsd))), genefilter::rowVars(assay(vsd)), main="VST")
}
\author{
  Simon Anders
}
\seealso{
  \code{\link{plotPCA}}, \code{\link{rlogTransformation}}
}

