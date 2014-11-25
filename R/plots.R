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
#' same manner and call \code{plotMA}.
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
#' @return A \code{trellis} object.
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
#' @param intgroup interesting groups: a character vector of names in \code{colData(x)} to use for grouping
#' @param ntop number of top genes to use for principal components, selected by highest
#'    row variance
#' @param returnData should the function only return the data.frame of PC1 and PC2
#' with intgroup covariates for custom plotting (default is FALSE)
#' 
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#' 
#' @author Wolfgang Huber
#'
#' @note See the vignette for an example of variance stabilization and PCA plots.
#' Note that the source code of plotPCA is very simple and commented.
#' Users should find it easy to customize this function.
#' 
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(betaSD=1)
#' rld <- rlog(dds)
#' plotPCA(rld)
#' 
#' @export
plotPCA = function(x, intgroup="condition", ntop=500, returnData=FALSE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(x))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(x)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  # add the intgroup factors together to create a new grouping factor
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop=FALSE])
  group <- factor(apply( intgroup.df, 1, paste, collapse=" : "))

  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(x))

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance"))
}


#' Plot of normalized counts for a single gene on log scale
#'
#' Note: normalized counts plus a pseudocount of 0.5 are shown.
#' 
#' @param dds a \code{DESeqDataSet}
#' @param gene a character, specifying the name of the gene to plot
#' @param intgroup interesting groups: a character vector of names in \code{colData(x)} to use for grouping
#' @param normalized whether the counts should be normalized by size factor
#' (default is TRUE)
#' @param transform whether to present log2 counts (TRUE) or
#' to present the counts on the log scale (FALSE, default)
#' @param main as in 'plot'
#' @param xlab as in 'plot'
#' @param returnData should the function only return the data.frame of counts and
#' covariates for custom plotting (default is FALSE)
#' @param ... arguments passed to plot
#' 
#' @examples
#'
#' dds <- makeExampleDESeqDataSet()
#' plotCounts(dds, "gene1")
#' 
#' @export
plotCounts <- function(dds, gene, intgroup="condition",
                       normalized=TRUE, transform=FALSE,
                       main, xlab="group",
                       returnData=FALSE, ...) {
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) & (gene >= 1 & gene <= nrow(dds)))))
  if (!all(intgroup %in% names(colData(dds)))) stop("all variables in 'intgroup' must be columns of colData")
  stopifnot(returnData | all(sapply(intgroup, function(v) is(colData(dds)[[v]], "factor"))))
  if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
    dds <- estimateSizeFactors(dds)
  }
  cnts <- counts(dds,normalized=normalized)[gene,]
  group <- if (length(intgroup) == 1) {
    colData(dds)[[intgroup]]
  } else if (length(intgroup) == 2) {
    lvls <- as.vector(t(outer(levels(colData(dds)[[intgroup[1]]]),
                              levels(colData(dds)[[intgroup[2]]]),
                              function(x,y) paste(x,y,sep=" : "))))
    droplevels(factor(apply( as.data.frame(colData(dds)[, intgroup, drop=FALSE]),
                            1, paste, collapse=" : "), levels=lvls))
  } else {
    factor(apply( as.data.frame(colData(dds)[, intgroup, drop=FALSE]),
                 1, paste, collapse=" : "))
  }
  data <- data.frame(count=cnts + .5, group=as.integer(group))
  if (transform) {
    data$count <- log2(data$count)
    ylab <- expression(log[2]~count)
    logxy <- ""
  } else {
    ylab <- ifelse(normalized,"normalized count","count")
    logxy <- "y"
  }
  if (missing(main)) {
    main <- if (is.numeric(gene)) {
      rownames(dds)[gene]
    } else {
      gene
    }
  }
  if (returnData) return(data.frame(count=data$count, colData(dds)[intgroup]))
  plot(data$group + runif(ncol(dds),-.05,.05), data$count, xlim=c(.5,max(data$group)+.5),
       log=logxy, xaxt="n", xlab=xlab, ylab=ylab, main=main, ...)
  axis(1, at=seq_along(levels(group)), levels(group))
}


##############
# unexported #
##############


# convenience function for adding alpha transparency to named colors
col2useful <- function(col,alpha) {
  x <- col2rgb(col)/255
  rgb(x[1],x[2],x[3],alpha)
}
