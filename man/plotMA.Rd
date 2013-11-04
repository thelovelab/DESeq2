\docType{methods}
\name{plotMA}
\alias{plotMA}
\alias{plotMA,DESeqDataSet-method}
\title{MA-plot from base means and log fold changes}
\usage{
  \S4method{plotMA}{DESeqDataSet}(object, alpha, main,
    ylim, ...)
}
\arguments{
  \item{object}{a DESeqDataSet processed by
  \code{\link{DESeq}}, or the individual functions
  \code{\link{nbinomWaldTest}} or \code{\link{nbinomLRT}}}

  \item{alpha}{the significance level for thresholding
  adjusted p-values}

  \item{main}{optional title for the plot}

  \item{ylim}{optional y limits}

  \item{...}{further arguments passed to
  \code{\link{results}}}
}
\description{
  A simple helper function that makes a so-called
  "MA-plot", i.e. a scatter plot of log2 fold changes (on
  the y-axis) versus the mean of normalized counts (on the
  x-axis). Note: it takes a \code{DESeqDataSet} object as
  its first argument.
}
\details{
  The function is practically 2 lines of code: building a
  dataframe and passing this to the plotMA method for
  \code{data.frame} from the geneplotter package. The code
  of this function can be seen with:
  \code{getMethod("plotMA","DESeqDataSet")} If users wish
  to modify the graphical parameters of the plot, it is
  recommended to build the dataframe in the same manner.
}
\examples{
dds <- makeExampleDESeqDataSet()
dds <- DESeq(dds)
plotMA(dds)
}
\author{
  Michael Love
}

