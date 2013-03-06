\docType{methods}
\name{normalizationFactors}
\alias{normalizationFactors}
\alias{normalizationFactors,DESeqSummarizedExperiment-method}
\alias{normalizationFactors<-}
\alias{normalizationFactors<-,DESeqSummarizedExperiment,matrix-method}
\title{Accessor functions for the normalization factors in a DESeqSummarizedExperiment
object.}
\usage{
  \S4method{normalizationFactors}{DESeqSummarizedExperiment}(object)

  \S4method{normalizationFactors}{DESeqSummarizedExperiment,matrix}(object)<-value
}
\arguments{
  \item{object}{a \code{DESeqSummarizedExperiment} object.}

  \item{value}{the matrix of normalization factors}
}
\description{
  Gene-specific normalization factors for each sample can
  be provided as a matrix, which will preempt
  \code{\link{sizeFactors}}. In some experiments, counts
  for each sample have varying dependence on covariates,
  e.g. on GC-content for sequencing data run on different
  days, and in this case it makes sense to provide
  gene-specific factors for each sample rather than a
  single size factor.
}
\details{
  Normalization factors alter the model of
  \code{\link{DESeq}} in the following way, for counts
  \eqn{K_{ij}}{K_ij} for gene i and sample j:

  \deqn{K_{ij} \sim \textrm{NB}(NF_{ij} \mu_{ij},
  \alpha_i)}{K_ij ~ NB(NF_ij * mu_ij, alpha_i)}
}
\note{
  Normalization factors are on the scale of the counts
  (similar to \code{\link{sizeFactors}}) and unlike
  offsets, which are typically on the scale of the
  predictors (in this case, log counts). Normalization
  factors should include size factor normalization and
  should have a mean around 1, as is the case with size
  factors.
}
\examples{
dse <- makeExampleDESeqSummarizedExperiment()
normFactors <- matrix(runif(nrow(dse)*ncol(dse),0.5,1.5),
                      ncol=ncol(dse),nrow=nrow(dse))
normalizationFactors(dse) <- normFactors
dse <- estimateDispersions(dse)
dse <- nbinomWaldTest(dse)
}

