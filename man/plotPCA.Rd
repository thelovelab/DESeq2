\name{plotPCA}
\alias{plotPCA}
\title{Sample PCA plot from variance-stabilized data}
\usage{
  plotPCA(x, intgroup = "condition", ntop = 500)
}
\arguments{
  \item{x}{a SummarizedExperiment, with transformed data in
  \code{assay(x)}, produced by
  \code{\link{varianceStabilizingTransformation}}}

  \item{intgroup}{a character vector of names in
  \code{colData(x)} to use for grouping}

  \item{ntop}{number of top genes to use for principal
  components, selected by highest row variance}
}
\description{
  This plot helps to check for batch effects and the like.
}
\note{
  See the vignette for an example of variance stabilization
  and PCA plots.
}
\examples{
dds <- makeExampleDESeqDataSet(betaSd=1)
design(dds) <- formula(~ 1)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd)
}
\author{
  Wolfgang Huber
}

