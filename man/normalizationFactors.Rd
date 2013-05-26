\docType{methods}
\name{normalizationFactors}
\alias{normalizationFactors}
\alias{normalizationFactors,DESeqDataSet-method}
\alias{normalizationFactors<-}
\alias{normalizationFactors<-,DESeqDataSet,matrix-method}
\title{Accessor functions for the normalization factors in a DESeqDataSet
object.}
\usage{
  \S4method{normalizationFactors}{DESeqDataSet}(object)

  \S4method{normalizationFactors}{DESeqDataSet,matrix}(object)<-value
}
\arguments{
  \item{object}{a \code{DESeqDataSet} object.}

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
  \eqn{K_{ij}}{K_ij} and normalization factors
  \eqn{NF_{ij}}{NF_ij} for gene i and sample j:

  \deqn{ K_{ij} \sim \textrm{NB}( \mu_{ij}, \alpha_i) }{
  K_ij ~ NB(mu_ij, alpha_i) } \deqn{ \mu_{ij} = NF_{ij}
  q_{ij} }{ mu_ij = NF_ij * q_ij }
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
dds <- makeExampleDESeqDataSet()
normFactors <- matrix(runif(nrow(dds)*ncol(dds),0.5,1.5),
                      ncol=ncol(dds),nrow=nrow(dds))
normalizationFactors(dds) <- normFactors
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
}

