\name{plotMA}
\alias{plotMA}
\title{MA-plot from base means and log fold changes}
\usage{
  plotMA(dds, lfcColname, pvalColname, pvalCutoff = 0.1,
    ylim,
    col = ifelse(mcols(dds)[, pvalColname] < pvalCutoff, "red", "black"),
    linecol = "#ff000080",
    xlab = "mean of normalized counts",
    ylab = expression(log[2] ~ fold ~ change), log = "x",
    cex = 0.45, ...)
}
\arguments{
  \item{dds}{a DESeqDataSet}

  \item{lfcColname}{the name of the column for log fold
  changes, if not provided this will default to the last
  variable in the design formula}

  \item{pvalColname}{the name of the column for
  pvalues/adjusted pvalues, if not provided this will
  default to \code{WaldAdjPvalue_lfcColname}}

  \item{pvalCutoff}{the cutoff for drawing red or black
  points}

  \item{ylim}{ylim}

  \item{col}{col}

  \item{linecol}{the color of the horizontal line}

  \item{xlab}{xlab}

  \item{ylab}{ylab}

  \item{log}{log, defaults to "x", the y-axis is already in
  log scale}

  \item{cex}{cex}

  \item{...}{further arguments to \code{plot}}
}
\description{
  A simple helper function that makes a so-called
  "MA-plot", i.e. a scatter plot of log2 fold changes (on
  the y-axis) versus the mean of normalized counts (on the
  x-axis).
}
\examples{
dds <- makeExampleDESeqDataSet()
dds <- DESeq(dds)
plotMA(dds)
}
\author{
  Wolfgang Huber
}

