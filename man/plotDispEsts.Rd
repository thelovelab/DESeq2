\name{plotDispEsts}
\alias{plotDispEsts}
\title{Plot dispersion estimates}
\usage{
  plotDispEsts(dds, ymin, genecol = "black",
    fitcol = "red", finalcol = "dodgerblue", legend = TRUE,
    xlab = "mean of normalized counts",
    ylab = "dispersion", log = "xy", cex = 0.45, ...)
}
\arguments{
  \item{dds}{a DESeqDataSet}

  \item{ymin}{the lower bound for points on the plot,
  points beyond this are drawn as triangles at ymin}

  \item{genecol}{the color for gene-wise dispersion
  estimates}

  \item{fitcol}{the color of the fitted estimates}

  \item{finalcol}{the color of the final estimates used for
  testing}

  \item{legend}{logical, whether to draw a legend}

  \item{xlab}{xlab}

  \item{ylab}{ylab}

  \item{log}{log}

  \item{cex}{cex}

  \item{...}{further arguments to \code{plot}}
}
\description{
  A simple helper function that plots the per-gene
  dispersion estimates together with the fitted
  mean-dispersion relationship.
}
\examples{
dds <- makeExampleDESeqDataSet()
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
plotDispEsts(dds)
}
\author{
  Simon Anders
}

