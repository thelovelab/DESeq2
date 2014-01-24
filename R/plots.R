#' Plot dispersion estimates
#'
#' A simple helper function that plots the per-gene dispersion
#' estimates together with the fitted mean-dispersion relationship.
#'
#' @usage
#' \S4method{plotDispEsts}{DESeqDataSet}(object, ymin,
#'   genecol = "black", fitcol = "red", finalcol = "dodgerblue",
#'   legend=TRUE, xlab, ylab, log = "xy", cex = 0.45, ...)
#'
#' @docType methods
#' @name plotDispEsts
#' @rdname plotDispEsts
#' @aliases plotDispEsts plotDispEsts,DESeqDataSet-method
#' 
#' @param object a DESeqDataSet
#' @param ymin the lower bound for points on the plot, points beyond this
#'    are drawn as triangles at ymin
#' @param genecol the color for gene-wise dispersion estimates
#' @param fitcol the color of the fitted estimates
#' @param finalcol the color of the final estimates used for testing
#' @param legend logical, whether to draw a legend
#' @param xlab xlab
#' @param ylab ylab
#' @param log log
#' @param cex cex
#' @param ... further arguments to \code{plot}
#'
#' @author Simon Anders
#'
#' @examples
#' 
#' dds <- makeExampleDESeqDataSet()
#' dds <- estimateSizeFactors(dds)
#' dds <- estimateDispersions(dds)
#' plotDispEsts(dds)
plotDispEsts.DESeqDataSet <- function( object, ymin,
  genecol = "black", fitcol = "red", finalcol = "dodgerblue",
  legend=TRUE, xlab, ylab, log = "xy", cex = 0.45, ... )
{
  if (missing(xlab)) xlab <- "mean of normalized counts"
  if (missing(ylab)) ylab <- "dispersion"
  
  px = mcols(object)$baseMean
  sel = (px>0)
  px = px[sel]

  py = mcols(object)$dispGeneEst[sel]
  if(missing(ymin))
      ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)

  plot(px, pmax(py, ymin), xlab=xlab, ylab=ylab,
    log=log, pch=ifelse(py<ymin, 6, 20), col=col2useful(genecol,.8), cex=cex, ... )

  # use a circle over outliers
  pchOutlier <- ifelse(mcols(object)$dispOutlier[sel],1,16)
  cexOutlier <- ifelse(mcols(object)$dispOutlier[sel],2*cex,cex)
  lwdOutlier <- ifelse(mcols(object)$dispOutlier[sel],2,1)
  if (!is.null(dispersions(object))) {
    points(px, dispersions(object)[sel], col=col2useful(finalcol,.8), cex=cexOutlier,
           pch=pchOutlier, lwd=lwdOutlier)
  }

  if (!is.null(mcols(object)$dispFit)) {
    points(px, mcols(object)$dispFit[sel], col=col2useful(fitcol,.8), cex=cex, pch=16)
  }
  
  if (legend) {
    legend("bottomright",c("gene-est","fitted","final"),pch=16,
           col=c(genecol,fitcol,finalcol),bg="white")
  }
}

#' @rdname plotDispEsts
#' @export
setMethod("plotDispEsts", signature(object="DESeqDataSet"), plotDispEsts.DESeqDataSet)


#' MA-plot from base means and log fold changes
#'
#' A simple helper function that makes a so-called "MA-plot", i.e. a
#' scatter plot of log2 fold changes (on the y-axis) versus the mean of
#' normalized counts (on the x-axis).
#'
#' This function is essentially two lines of code: building a
#' \code{data.frame} and passing this to the \code{plotMA} method
#' for \code{data.frame} from the geneplotter package.
#' The code of this function can be seen with:
#' \code{getMethod("plotMA","DESeqDataSet")}
#' If users wish to modify the graphical parameters of the plot,
#' it is recommended to build the data.frame in the
#' same manner and call \code{plotMA}
#'
#' @usage
#' \S4method{plotMA}{DESeqResults}(object, alpha, main, ylim, ...)
#' \S4method{plotMA}{DESeqDataSet}(object, alpha, main, ylim, ...)
#' 
#'
#' @docType methods
#' @name plotMA
#' @rdname plotMA
#' @aliases plotMA plotMA,DESeqDataSet-method plotMA,DESeqResults-method
#' 
#' @param object a \code{DESeqResults} object produced by \code{\link{results}};
#' or a \code{DESeqDataSet} processed by \code{\link{DESeq}}, or the
#' individual functions \code{\link{nbinomWaldTest}} or \code{\link{nbinomLRT}}
#' @param alpha the significance level for thresholding adjusted p-values
#' @param main optional title for the plot
#' @param ylim optional y limits
#' @param ... further arguments passed to \code{plotMA} if object
#' is \code{DESeqResults} or to \code{\link{results}} if object is
#' \code{DESeqDataSet}
#'
#' @author Michael Love
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet()
#' dds <- DESeq(dds)
#' plotMA(dds)
#' res <- results(dds)
#' plotMA(res)
#'
#' @importFrom geneplotter plotMA
plotMA.DESeqDataSet <- function(object, alpha=.1, main="", ylim, ...) {
    res <- results(object, ...)
    plotMA.DESeqResults(res, alpha=alpha, main=main, ylim=ylim)
}

#' @rdname plotMA
#' @export
setMethod("plotMA", signature(object="DESeqDataSet"), plotMA.DESeqDataSet)

plotMA.DESeqResults <- function(object, alpha=.1, main="", ylim, ...) {
    df <- data.frame(mean = object$baseMean,
                     lfc = object$log2FoldChange,
                     isDE = ifelse(is.na(object$padj), FALSE, object$padj < alpha))
    if (missing(ylim)) {
      plotMA(df, main=main, ...)
    } else {
      plotMA(df, main=main, ylim=ylim, ...)
    }  
}

#' @rdname plotMA
#' @export
setMethod("plotMA", signature(object="DESeqResults"), plotMA.DESeqResults)


#' Sample PCA plot from variance-stabilized data
#' 
#' This plot helps to check for batch effects and the like.
#' 
#' @param x a SummarizedExperiment, with data in \code{assay(x)},
#' produced for example by either \code{\link{varianceStabilizingTransformation}}
#' or \code{\link{rlogTransformation}}
#' @param intgroup a character vector of names in \code{colData(x)} to use for grouping
#' @param ntop number of top genes to use for principal components, selected by highest
#'    row variance
#'
#' @author Wolfgang Huber
#'
#' @note See the vignette for an example of variance stabilization and PCA plots.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom lattice xyplot draw.key
#' 
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(betaSD=1)
#' vsd <- varianceStabilizingTransformation(dds)
#' plotPCA(vsd)
#'
#' @export
plotPCA = function(x, intgroup="condition", ntop=500)
{
  rv = rowVars(assay(x))
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(x)[select,]))

  fac = factor(apply( as.data.frame(colData(x)[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  colours = if( nlevels(fac) >= 3 )
    brewer.pal(nlevels(fac), "Paired")
  else  
    c( "lightgreen", "dodgerblue" )

  xyplot(PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=2,
    aspect = "iso", col=colours,
    main = draw.key(key = list(
      rect = list(col = colours),
      text = list(levels(fac)),
      rep = FALSE)))
}


# convenience function for adding alpha transparency to named colors
col2useful <- function(col,alpha) {
  x <- col2rgb(col)/255
  rgb(x[1],x[2],x[3],alpha)
}
