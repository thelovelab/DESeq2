\name{rlogTransformation}
\alias{rlogData}
\alias{rlogTransformation}
\title{Apply a 'regularized log' transformation}
\usage{
  rlogTransformation(object, unsupervised = TRUE,
    samplesVector, priorSigmasq, rowVarQuantile = 0.9)

  rlogData(object, samplesVector, priorSigmasq,
    rowVarQuantile = 0.9)
}
\arguments{
  \item{object}{a DESeqDataSet}

  \item{unsupervised}{logical, whether the dispersions
  should be re-estimated using a design formula with only
  the intercept. This is recommended in order to ensure
  that the data transformation is an unsupervised method.}

  \item{samplesVector}{a character vector or factor of the
  sample identifiers}

  \item{priorSigmasq}{a single value, the variance of the
  prior on the sample betas, which if missing is estimated
  from the rows which do not have any zeros}

  \item{rowVarQuantile}{the quantile of the row variances
  of log fold changes which will be used to set the width
  of the prior}
}
\value{
  for \code{rlogTransformation}, a SummarizedExperiment
  with assay data elements equal to \eqn{\log_2(q_{ij}) =
  X_{j.} \beta_i}{log2(q_ij) = X_j. * beta_i}, see formula
  at \code{\link{DESeq}}. for \code{rlogData}, a
  \code{matrix} of the same dimension as the count data,
  containing the transformed values.
}
\description{
  This function uses Tikhonov/ridge regularization, as in
  \code{\link{nbinomWaldTest}}, to transform the data to
  the log2 scale in a way which minimizes differences
  between samples for rows with small counts. The
  transformation produces a similar variance stabilizing
  effect as
  \code{\link{varianceStabilizingTransformation}}, though
  \code{rlogTransformation} is more robust in the case when
  the size factors vary widely. The transformation is
  useful when checking for outliers or as input for machine
  learning techniques such as clustering or linear
  discriminant analysis.
}
\details{
  The 'regularization' referred to here corresponds to the
  maximum a posteriori solution to the GLM with a prior on
  the coefficients for each sample. The prior width is
  calculated as follows: coefficients are fit for a model
  with a term for each sample and for the intercept. This
  would typically result in an unidentifiable solution, so
  a very wide prior is used. Then the prior variance is
  estimated by taking the mean of the row-wise variance of
  the sample coefficients. A second and final GLM fit is
  performed using this prior. It is also possible to supply
  the variance of the prior. See the vignette for an
  example of the use and a comparison with
  \code{varianceStabilizingTransformation}
}
\examples{
dds <- makeExampleDESeqDataSet(betaSd=1)
rld <- rlogTransformation(dds)
dists <- dist(t(assay(rld)))
plot(hclust(dists))
}
\seealso{
  \code{\link{plotPCA}},
  \code{\link{varianceStabilizingTransformation}}
}

