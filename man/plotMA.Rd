\docType{methods}
\name{plotMA}
\alias{plotMA}
\alias{plotMA,DESeqDataSet-method}
\title{MA-plot from base means and log fold changes}
\usage{
  \S4method{plotMA}{DESeqDataSet}(object, lfcColname,
    pvalues, pvalCutoff=.1, ylim, linecol = "#ff000080",
    pointcol = c("black","red"), xlab, ylab, log = "x",
    cex=0.45, ...)
}
\arguments{
  \item{object}{a DESeqDataSet processed by
  \code{\link{DESeq}}, or the individual functions
  \code{\link{nbinomWaldTest}} or \code{\link{nbinomLRT}}}

  \item{lfcColname}{the name of the column for log fold
  changes, if not provided this will default to the last
  variable in the design formula.  for options for this
  argument, check resultsNames(object).}

  \item{pvalues}{a vector of the p-values or adjusted
  p-values to use in coloring the points. If not provided,
  defaults to the 'padj' column of results(object)}

  \item{pvalCutoff}{the cutoff for drawing red or black
  points}

  \item{ylim}{the limits for the y axis (chosen
  automatically if not specified)}

  \item{linecol}{the color of the horizontal line}

  \item{pointcol}{a vector of length two of the colors for
  the not significant and significant points, respectively}

  \item{xlab}{the label for the x axis}

  \item{ylab}{the label for the y axis}

  \item{log}{log, defaults to "x", the y-axis is already in
  log scale}

  \item{cex}{the size of the points}

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

