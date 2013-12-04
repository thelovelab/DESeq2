\name{rlogTransformation}
\alias{rlogData}
\alias{rlogTransformation}
\title{Apply a 'regularized log' transformation}
\usage{
  rlogTransformation(object, blind = TRUE, samplesVector,
    betaPriorVar, intercept)

  rlogData(object, samplesVector, betaPriorVar, intercept)
}
\arguments{
  \item{object}{a DESeqDataSet}

  \item{blind}{logical, whether to blind the transformation
  to the experimental design. blind=TRUE should be used for
  comparing samples in an manner unbiased by prior
  information on samples, for example to perform sample QA
  (quality assurance). blind=FALSE should be used for
  transforming data for downstream analysis, where the full
  use of the design information should be made.}

  \item{samplesVector}{a character vector or factor of the
  sample identifiers}

  \item{betaPriorVar}{a single value, the variance of the
  prior on the sample betas, which if missing is estimated
  from the data}

  \item{intercept}{by default, this is not provided and
  calculated automatically. if provided, this should be a
  vector as long as the number of rows of object, which is
  log2 of the mean normalized counts from a previous
  dataset. this will enforce the intercept for the GLM,
  allowing for a "frozen" rlog transformation based on a
  previous dataset.}
}
\value{
  for \code{rlogTransformation}, a SummarizedExperiment
  with assay data elements equal to \eqn{\log_2(q_{ij}) =
  x_{j.} \beta_i}{log2(q_ij) = x_j. * beta_i}, see formula
  at \code{\link{DESeq}}. for \code{rlogData}, a
  \code{matrix} of the same dimension as the count data,
  containing the transformed values. To avoid returning
  matrices with NA values where there were zeros for all
  rows of the unnormalized counts, rlogTransformation
  returns instead all zeros (essentially adding a
  pseudocount of one, only to those rows in which all
  samples have zero).
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
  the coefficients for each sample. The fitted dispersions
  are used rather than the MAP dispersions (so similar to
  the \code{\link{varianceStabilizingTransformation}}) as
  the blind dispersion estimation would otherwise shrink
  large, true log fold changes. The prior variance is
  calculated as follows: a matrix is constructed of the
  logarithm of the counts plus a pseudocount of 0.5, the
  row means of these log counts are then subtracted,
  leaving an estimate of the log fold changes per sample.
  The prior variance is set to the variance of all log fold
  change estimates. A second and final GLM fit is then
  performed using this prior. It is also possible to supply
  the variance of the prior. See the vignette for an
  example of the use and a comparison with
  \code{varianceStabilizingTransformation}

  The parameters of the rlog transformation from a previous
  dataset can be "frozen" and reapplied to new samples. See
  the "Data quality assessment" section of the vignette for
  strategies to see if new samples are sufficiently similar
  to previous datasets. The "freezing" is accomplished by
  saving the dispersion function, beta prior variance and
  the intercept from a previous dataset, and running
  rlogTransformation with 'blind' set to FALSE (see example
  below).
}
\examples{
dds <- makeExampleDESeqDataSet(betaSD=1)
rld <- rlogTransformation(dds, blind=TRUE)
dists <- dist(t(assay(rld)))
plot(hclust(dists))

# run the rlog transformation on one dataset
design(dds) <- ~ 1
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
rld <- rlogTransformation(dds, blind=FALSE)

# apply the parameters to a new sample

ddsNew <- makeExampleDESeqDataSet(m=1)
mcols(ddsNew)$dispFit <- mcols(dds)$dispFit
betaPriorVar <- attr(rld,"betaPriorVar")
intercept <- mcols(rld)$rlogIntercept
rldNew <- rlogTransformation(ddsNew, blind=FALSE,
                           betaPriorVar=betaPriorVar,
                           intercept=intercept)
}
\seealso{
  \code{\link{plotPCA}},
  \code{\link{varianceStabilizingTransformation}}
}

