# convenience function for adding alpha transparency to named colors
col2useful <- function(col,alpha) {
  x <- col2rgb(col)/255
  rgb(x[1],x[2],x[3],alpha)
}


#' Plot dispersion estimates
#'
#' A simple helper function that plots the per-gene dispersion
#' estimates together with the fitted mean-dispersion relationship.
#'
#' @param dse a DESeqSummarizedExperiment
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
#' dse <- makeExampleDESeqSummarizedExperiment()
#' dse <- estimateSizeFactors(dse)
#' dse <- estimateDispersions(dse)
#' plotDispEsts(dse)
#' 
#' @export
plotDispEsts = function( dse, ymin,
  genecol = "black", fitcol = "red", finalcol = "blue",
  legend=TRUE, xlab = "mean of normalized counts",
  ylab = "dispersion", log = "xy", cex = 0.45, ... )
{ 
  px = mcols(dse)$baseMean
  sel = (px>0)
  px = px[sel]

  py = mcols(dse)$dispGeneEst[sel]
  if(missing(ymin))
      ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)

  plot(px, pmax(py, ymin), xlab=xlab, ylab=ylab,
    log=log, pch=ifelse(py<ymin, 6, 20), col=genecol, cex=cex, ... )
  points(px, dispersions(dse)[sel], col=col2useful(finalcol,.3), cex=cex, pch=20)
  points(px, mcols(dse)$dispFit[sel], col=col2useful(fitcol,.5), cex=cex, pch=20)
  if (legend) {
    legend("bottomright",c("gene-est","fitted","final"),pch=20,col=c(genecol,fitcol,finalcol),bg="white")
  }
}


#' MA-plot from base means and log fold changes
#'
#' A simple helper function that makes a so-called "MA-plot", i.e. a
#' scatter plot of log2 fold changes (on the y-axis) versus the mean of
#' normalized counts (on the x-axis).
#'
#' @param dse a DESeqSummarizedExperiment
#' @param lfcColname the name of the column for log fold changes, if
#'    not provided this will default to the last variable in the design formula
#' @param pvalColname the name of the column for pvalues/adjusted pvalues, if
#'    not provided this will default to \code{WaldAdjPvalue_lfcColname}
#' @param pvalCutoff the cutoff for drawing red or black points
#' @param ylim ylim
#' @param col col
#' @param linecol the color of the horizontal line
#' @param xlab xlab
#' @param ylab ylab
#' @param log log, defaults to "x", the y-axis is already in log scale
#' @param cex cex
#' @param ... further arguments to \code{plot}
#'
#' @author Wolfgang Huber
#'
#' @examples
#'
#' dse <- makeExampleDESeqSummarizedExperiment()
#' dse <- DESeq(dse)
#' plotMA(dse)
#' 
#' @export
plotMA = function(dse, lfcColname, pvalColname, pvalCutoff=.1, ylim, col = ifelse(mcols(dse)[,pvalColname] < pvalCutoff, "red", "black"), linecol = "#ff000080", xlab = "mean of normalized counts", ylab = expression(log[2]~fold~change), log = "x", cex=0.45, ...)
{
  # if not specified, try the last variable of the design formula,
  # the last level of this variable, and the Wald test adjusted p-values
  if (missing(lfcColname)) {
    lfcColname <- lastCoefName(dse)
  }
  if (missing(pvalColname)) {
    pvalColname <- paste0("WaldAdjPvalue_",lfcColname)
  }
  x = mcols(dse)
  col = col[x$baseMean > 0]
  x = x[x$baseMean > 0,]
  py = x[,lfcColname]
  if(missing(ylim))
      ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
  plot(x$baseMean, pmax(ylim[1], pmin(ylim[2], py)),
       log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 20)),
       cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
  abline(h=0, lwd=4, col=linecol)
}



#' Sample PCA plot from variance-stabilized data
#' 
#' This plot helps to check for batch effects and the like.
#' 
#' @param x a SummarizedExperiment, with transformed data in \code{assay(x)},
#' produced by \code{\link{varianceStabilizingTransformation}}
#' @param intgroup a character vector of names in \code{colData(x)} to use for grouping
#' @param ntop number of top genes to use for principal components, selected by highest
#'    row variance
#'
#' @author Wolfgang Huber
#'
#' @note See the vignette for an example of variance stabilization and PCA plots.
#'
#' @examples
#'
#' dse <- makeExampleDESeqSummarizedExperiment(betaSd=1)
#' design(dse) <- formula(~ 1)
#' dse <- estimateSizeFactors(dse)
#' dse <- estimateDispersions(dse)
#' vsd <- varianceStabilizingTransformation(dse)
#' plotPCA(vsd)
#'
#' @export
plotPCA = function(x, intgroup="condition", ntop=500)
{
  rv = rowVars(assay(x))
  select = order(rv, decreasing=TRUE)[seq_len(ntop)]
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
