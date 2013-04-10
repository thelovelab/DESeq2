\name{varianceStabilizingTransformation}
\alias{getVarianceStabilizedData}
\alias{varianceStabilizingTransformation}
\title{Apply a variance stabilizing transformation (VST) to the count data}
\usage{
  varianceStabilizingTransformation(object)

  getVarianceStabilizedData(object)
}
\arguments{
  \item{object}{a DESeqDataSet, with \code{design(object)
  <- formula(~ 1)} and size factors (or normalization
  factors) and dispersions estimated using local or
  parametric \code{fitType}.}
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
  adding the shot noise). The function requires a blind
  estimate of the variance function, i.e., one ignoring
  conditions. This is achieved by setting the design
  formula to \code{~ 1} and then calling
  \code{\link{estimateDispersions}}.

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
design(dds) <- formula(~ 1)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
vsd <- varianceStabilizingTransformation(dds)
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

