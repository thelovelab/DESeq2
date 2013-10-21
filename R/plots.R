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
#' 
#' @export
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
#' @usage
#' \S4method{plotMA}{DESeqDataSet}(object, lfcColname, pvalues, pvalCutoff=.1, ylim,
#'   linecol = "#ff000080", pointcol = c("black","red"),
#'   xlab, ylab, log = "x", cex=0.45, ...)
#'
#' @docType methods
#' @name plotMA
#' @rdname plotMA
#' @aliases plotMA plotMA,DESeqDataSet-method
#' 
#' @param object a DESeqDataSet processed by \code{\link{DESeq}}, or the
#' individual functions \code{\link{nbinomWaldTest}} or \code{\link{nbinomLRT}}
#' @param lfcColname the name of the column for log fold changes, if
#'    not provided this will default to the last variable in the design formula.
#'    for options for this argument, check resultsNames(object).
#' @param pvalues a vector of the p-values or adjusted p-values to use in coloring
#'    the points. If not provided, defaults to the 'padj' column of results(object)
#' @param pvalCutoff the cutoff for drawing red or black points
#' @param ylim the limits for the y axis (chosen automatically if not specified)
#' @param linecol the color of the horizontal line
#' @param pointcol a vector of length two of the colors for the not significant and
#' significant points, respectively
#' @param xlab the label for the x axis
#' @param ylab the label for the y axis
#' @param log log, defaults to "x", the y-axis is already in log scale
#' @param cex the size of the points
#' @param ... further arguments to \code{plot}
#'
#' @author Wolfgang Huber
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet()
#' dds <- DESeq(dds)
#' plotMA(dds)
#' 
#' @export
plotMA.DESeqDataSet <- function(object, lfcColname, pvalues,
  pvalCutoff=.1, ylim, linecol = "#ff000080",
  pointcol = c("black","red"), xlab, ylab, log = "x",
  cex=0.45, ...) {
  
  if (missing(xlab)) xlab <- "mean of normalized counts"
  if (missing(ylab)) ylab <- expression(log[2]~fold~change)
  
  if (!missing(pvalues)) {
    if (length(pvalues) != nrow(object)) {
      stop("length of pvalues should be equal to the number of rows of object")
    }
  }
  stopifnot(length(pointcol) == 2)

  if (!"results" %in% mcols(mcols(object))$type) {
    stop("first run DESeq() in order to produce an MA-plot")
  }
  # if not specified, try the last variable of the design formula,
  # the last level of this variable, and the Wald test adjusted p-values
  if (missing(lfcColname)) {
    lfcColname <- lastCoefName(object)
  }
  if (length(lfcColname) != 1 | !is.character(lfcColname)) {
    stop("the argument 'lfcColname' should be a character vector of length 1")
  }
  if (missing(pvalues)) {
    res <- results(object,name=lfcColname)
    pvalues <- res$padj
  }
  x <- mcols(object)
  
  stopifnot( length(cex) == 1 )
  col <- ifelse(is.na(pvalues) | pvalues > pvalCutoff, pointcol[1], pointcol[2])
  
  col = col[x$baseMean > 0]
  x = x[x$baseMean > 0,]
  py = x[,lfcColname]
  if (missing(ylim))
    ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
  plot(x$baseMean, pmax(ylim[1], pmin(ylim[2], py)),
       log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 20)),
       cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
  abline(h=0, lwd=4, col=linecol)
}

#' @rdname plotMA
#' @export
setMethod("plotMA", signature(object="DESeqDataSet"), plotMA.DESeqDataSet)

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
#' design(dds) <- formula(~ 1)
#' dds <- estimateSizeFactors(dds)
#' dds <- estimateDispersions(dds)
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
  if( nlevels(fac) >= 3 )
     colours = brewer.pal(nlevels(fac), "Paired")
  else  
   colours = c( "lightgreen", "dodgerblue" )

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
