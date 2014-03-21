\name{rlog}
\alias{rlog}
\alias{rlogData}
\alias{rlogDataFast}
\alias{rlogTransformation}
\title{Apply a 'regularized log' transformation}
\usage{
rlog(object, blind = TRUE, fast = FALSE, intercept, betaPriorVar, B)

rlogTransformation(object, blind = TRUE, fast = FALSE, intercept,
  betaPriorVar, B)

rlogData(object, intercept, betaPriorVar)

rlogDataFast(object, intercept, betaPriorVar, B)
}
\arguments{
  \item{object}{a DESeqDataSet}

  \item{blind}{logical, whether to blind the transformation
  to the experimental design. blind=TRUE should be used for
  comparing samples in an manner unbiased by prior
  information on samples, for example to perform sample QA
  (quality assurance). blind=FALSE should be used for
  transforming data for downstream analysis, where the full
  use of the design information should be made. If many of
  genes have large differences in counts due to the
  experimental design, it is important to set blind=FALSE
  for downstream analysis.}

  \item{fast}{if TRUE, an approximation to the rlog
  transformation is made, wherein an optimal amount of
  shrinkage of sample-wise log fold changes is calculated
  for genes divided into 12 batched based on mean counts.
  then shrinkage is performed by interpolating between
  these optimal values. optimal is defined by maximizing
  the same posterior as in the standard rlog
  transformation.}

  \item{intercept}{by default, this is not provided and
  calculated automatically. if provided, this should be a
  vector as long as the number of rows of object, which is
  log2 of the mean normalized counts from a previous
  dataset. this will enforce the intercept for the GLM,
  allowing for a "frozen" rlog transformation based on a
  previous dataset.}

  \item{betaPriorVar}{a single value, the variance of the
  prior on the sample betas, which if missing is estimated
  from the data}

  \item{B}{for the \code{fast} method only. by default,
  this is not provided and calculated automatically. if
  provided, this should be a vector as long as the number
  of rows of object. this is the amount of shrinkage from 0
  to 1 for each row, and takes precedence over internal
  calculation of B using betaPriorVar.}
}
\value{
\code{rlog}, equivalent to \code{rlogTransformation},
returns a \code{\link{SummarizedExperiment}}. The matrix of
transformed values are accessed by \code{assay(rld)}. for
\code{rlogData}, a \code{matrix} of the same dimension as
the count data, containing the transformed values. To avoid
returning matrices with NA values where there were zeros
for all rows of the unnormalized counts, rlogTransformation
returns instead all zeros (essentially adding a pseudocount
of one, only to those rows in which all samples have zero).
}
\description{
This function uses Tikhonov/ridge regularization to
transform the data to the log2 scale in a way which
minimizes differences between samples for rows with small
counts. The rlog transformation produces a similar variance
stabilizing effect as
\code{\link{varianceStabilizingTransformation}}, though
\code{rlogTransformation} is more robust in the case when
the size factors vary widely. The transformation is useful
when checking for outliers or as input for machine learning
techniques such as clustering or linear discriminant
analysis. Note: \code{rlog} is equivalent to
\code{rlogTransformation}, which take as input a
\code{\link{DESeqDataSet}} and return a
\code{\link{SummarizedExperiment}} object.
}
\details{
Note that neither rlog transformation nor the VST are used
by the differential expression estimation in
\code{\link{DESeq}}, which always occurs on the raw count
data, through generalized linear modeling which
incorporates knowledge of the variance-mean dependence. The
rlog transformation and VST are offered as separate
functionality which can be used for visualization,
clustering or other machine learning tasks. See the
transformation section of the vignette for more details.

The 'regularization' referred to here corresponds to the
maximum a posteriori solution to the GLM with a prior on
the coefficients for each sample. The fitted dispersions
are used rather than the MAP dispersions (so similar to the
\code{\link{varianceStabilizingTransformation}}) as the
blind dispersion estimation would otherwise shrink large,
true log fold changes. The prior variance is calculated as
follows: a matrix is constructed of the logarithm of the
counts plus a pseudocount of 0.5, the log of the row means
is then subtracted, leaving an estimate of the log fold
changes per sample over the fitted value using only an
intercept. The prior variance is then calculated as with
\code{\link{nbinomWaldTest}}, by matching the upper
quantiles of the observed log fold change estimates with an
upper quantile of the normal distribution. A second and
final GLM fit is then performed using this prior. It is
also possible to supply the variance of the prior. See the
vignette for an example of the use and a comparison with
\code{varianceStabilizingTransformation}.

The transformed values, rlog(K), are equal to
\eqn{rlog(K_{ij}) = \log_2(q_{ij}) = x_{j.}
\beta_i}{rlog(K_ij) = log2(q_ij) = x_j. * beta_i}, with
formula terms defined in \code{\link{DESeq}}.

The parameters of the rlog transformation from a previous
dataset can be "frozen" and reapplied to new samples. See
the "Data quality assessment" section of the vignette for
strategies to see if new samples are sufficiently similar
to previous datasets. The "freezing" is accomplished by
saving the dispersion function, beta prior variance and the
intercept from a previous dataset, and running
rlogTransformation with 'blind' set to FALSE (see example
below).
}
\examples{
dds <- makeExampleDESeqDataSet(m=6,betaSD=1)
rld <- rlog(dds, blind=TRUE)
dists <- dist(t(assay(rld)))
plot(hclust(dists))

# run the rlog transformation on one dataset
design(dds) <- ~ 1
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
rld <- rlog(dds, blind=FALSE)

# apply the parameters to a new sample

ddsNew <- makeExampleDESeqDataSet(m=1)
mcols(ddsNew)$dispFit <- mcols(dds)$dispFit
betaPriorVar <- attr(rld,"betaPriorVar")
intercept <- mcols(rld)$rlogIntercept
rldNew <- rlog(ddsNew, blind=FALSE,
               intercept=intercept,
               betaPriorVar=betaPriorVar)
}
\seealso{
\code{\link{plotPCA}},
\code{\link{varianceStabilizingTransformation}}
}

