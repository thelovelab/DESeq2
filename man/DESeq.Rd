\name{DESeq}
\alias{DESeq}
\title{Differential expression analysis based on the negative binomial distribution}
\usage{
  DESeq(object, fitType = c("parametric", "local", "mean"),
    betaPrior = TRUE, pAdjustMethod = "BH", quiet = FALSE)
}
\arguments{
  \item{object}{a DESeqDataSet object, see the constructor
  functions \code{\link{DESeqDataSet}},
  \code{\link{DESeqDataSetFromMatrix}},
  \code{\link{DESeqDataSetFromHTSeqCount}}.}

  \item{fitType}{either "parametric", "local", or "mean"
  for the type of fitting of dispersions to the mean
  intensity. See \code{\link{estimateDispersions}} for
  description.}

  \item{betaPrior}{whether or not to put a zero-mean normal
  prior on the non-intercept coefficients (Tikhonov/ridge
  regularization) See \code{\link{nbinomWaldTest}} for
  description.}

  \item{pAdjustMethod}{the method to use for adjusting
  p-values, see \code{?p.adjust}}

  \item{quiet}{whether to print messages at each step}
}
\value{
  a \code{\link{DESeqDataSet}} object with results stored
  as metadata columns.  The results can be accessed by
  calling the \code{\link{results}} function.  By default
  this will return the log2 fold changes and p-values for
  the last variable in the design formula.  See
  \code{\link{results}} for how to access results for other
  variables.
}
\description{
  This function performs a default analysis by calling, in
  order, the functions: \code{\link{estimateSizeFactors}},
  \code{\link{estimateDispersions}},
  \code{\link{nbinomWaldTest}}.
}
\details{
  The differential expression analysis uses a generalized
  linear model of the form:

  \deqn{ K_{ij} \sim \textrm{NB}( \mu_{ij}, \alpha_i) }{
  K_ij ~ NB(mu_ij, alpha_i) } \deqn{ \mu_{ij} = s_j q_{ij}
  }{ mu_ij = s_j * q_ij } \deqn{ \log_2(q_{ij}) = x_{j.}
  \beta_i }{ log2(q_ij) = x_j. * beta_i }

  where counts \eqn{K_{ij}}{K_ij} for gene i, sample j are
  modeled using a negative binomial distribution with
  fitted mean \eqn{\mu_{ij}}{mu_ij} and a gene-specific
  dispersion parameter \eqn{\alpha_i}{alpha_i}. The fitted
  mean is composed of a sample-specific size factor
  \eqn{s_j}{s_j} and a parameter \eqn{q_{ij}}{q_ij}
  proportional to the expected true concentration of
  fragments for sample j. The coefficients
  \eqn{\beta_i}{beta_i} give the log2 fold changes for gene
  i for each column of the model matrix \eqn{X}{X}.  The
  sample-specific size factors can be replaced by
  gene-specific normalization factors for each sample using
  \code{\link{normalizationFactors}}.  For details on the
  fitting of the log2 fold changes and calculation of
  p-values see \code{\link{nbinomWaldTest}}.
}
\examples{
dds <- makeExampleDESeqDataSet(betaSd=1)
dds <- DESeq(dds)
res <- results(dds)
}
\author{
  Michael Love
}
\references{
  Simon Anders, Wolfgang Huber: Differential expression
  analysis for sequence count data. Genome Biology 11
  (2010) R106,
  \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
}
\seealso{
  \code{\link{nbinomWaldTest}}, \code{\link{nbinomLRT}}
}

