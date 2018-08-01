#' Unmix samples using loss in a variance stabilized space
#'
#' Unmixes samples in \code{x} according to \code{pure} components,
#' using numerical optimization. The components in \code{pure}
#' are added on the scale of gene expression (either normalized counts, or TPMs).
#' The loss function when comparing fitted expression to the
#' samples in \code{x} occurs in a variance stabilized space.
#' This task is sometimes referred to as "deconvolution",
#' and can be used, for example, to identify contributions from
#' various tissues.
#' Note: some groups have found that the mixing contributions
#' may be more accurate if very lowly expressed genes across \code{x}
#' and \code{pure} are first removed. We have not explored this fully.
#' Note: if the \code{pbapply} package is installed a progress bar
#' will be displayed while mixing components are fit.
#'
#' @param x normalized counts or TPMs of the samples to be unmixed
#' @param pure normalized counts or TPMs of the "pure" samples
#' @param alpha for normalized counts, the dispersion of the data
#' when a negative binomial model is fit. this can be found by examining
#' the asymptotic value of \code{dispersionFunction(dds)}, when using
#' \code{fitType="parametric"} or the mean value when using
#' \code{fitType="mean"}.
#' @param shift for TPMs, the shift which approximately stabilizes the variance
#' of log shifted TPMs. Can be assessed with \code{vsn::meanSdPlot}.
#' @param power either 1 (for L1) or 2 (for squared) loss function.
#' Default is 1.
#' @param format \code{"matrix"} or \code{"list"}, default is \code{"matrix"}.
#' whether to output just the matrix of mixture components, or a list (see Value).
#' 
#' @param quiet suppress progress bar. default is FALSE, show progress bar
#' if pbapply is installed.
#'
#' @return a matrix, the mixture components for each sample in \code{x} (rows).
#' The "pure" samples make up the columns, and so each row sums to 1.
#' If colnames existed on the input matrices they will be propagated to the output matrix.
#' If \code{format="list"}, then a list, containing as elements:
#' (1) the matrix of mixture components,
#' (2) the correlations in the variance stabilized space of the fitted samples
#' to the samples in \code{x}, and
#' (3) the fitted samples as a matrix with the same dimension as \code{x}.
#'
#' @examples
#'
#' # some artificial data
#' cts <- matrix(c(80,50,1,100,
#'                 1,1,60,100,
#'                 0,50,60,100), ncol=4, byrow=TRUE)
#' # make a DESeqDataSet
#' dds <- DESeqDataSetFromMatrix(cts,
#'   data.frame(row.names=seq_len(ncol(cts))), ~1)
#' colnames(dds) <- paste0("sample",1:4)
#'
#' # note! here you would instead use
#' # estimateSizeFactors() to do actual normalization
#' sizeFactors(dds) <- rep(1, ncol(dds))
#'
#' norm.cts <- counts(dds, normalized=TRUE)
#'
#' # 'pure' should also have normalized counts...
#' pure <- matrix(c(10,0,0,
#'                  0,0,10,
#'                  0,10,0), ncol=3, byrow=TRUE)
#' colnames(pure) <- letters[1:3]
#' 
#' # for real data, you need to find alpha after fitting estimateDispersions()
#' mix <- unmix(norm.cts, pure, alpha=0.01)
#' 
#' @export
unmix <- function(x, pure, alpha, shift, power=1, format="matrix", quiet=FALSE) {

  format <- match.arg(format, c("matrix","list"))
  if (missing(alpha)) stopifnot(!missing(shift))
  if (missing(shift)) stopifnot(!missing(alpha))
  stopifnot(missing(shift) | missing(alpha))
  stopifnot(power %in% 1:2)
  stopifnot(nrow(x) == nrow(pure))
  stopifnot(ncol(pure) > 1)
  
  if (requireNamespace("pbapply", quietly=TRUE) & !quiet) {
    lapply <- pbapply::pblapply
  }

  cor.msg <- "some columns of 'pure' are highly correlated (>.99 after VST),
  may result in instabilty of unmix(). visually inspect correlations of 'pure'"
  
  if (missing(shift)) {
    stopifnot(alpha > 0)
    # variance stabilizing transformation for NB w/ fixed dispersion alpha
    vst <- function(q, alpha) ( 2 * asinh(sqrt(alpha * q)) - log(alpha) - log(4) ) / log(2)
    pure.cor <- cor(vst(pure, alpha)); diag(pure.cor) <- 0
    if (any(pure.cor > .99)) warning(cor.msg)
    sumLossVST <- function(p, i, vst, alpha, power) {
      sum(abs(vst(x[,i], alpha) - vst(pure %*% p, alpha))^power)
    }
    res <- lapply(seq_len(ncol(x)), function(i) {
      optim(par=rep(1, ncol(pure)), fn=sumLossVST, gr=NULL, i, vst, alpha, power,
            method="L-BFGS-B", lower=0, upper=100)$par
    })
  } else {
    stopifnot(shift > 0)
    # VST of shifted log
    vstSL <- function(q, shift) log(q + shift)
    pure.cor <- cor(vstSL(pure, shift)); diag(pure.cor) <- 0
    if (any(pure.cor > .99)) warning(cor.msg)
    sumLossSL <- function(p, i, vst, shift, power) {
      sum(abs(vstSL(x[,i], shift) - vstSL(pure %*% p, shift))^power)
    }
    res <- lapply(seq_len(ncol(x)), function(i) {
      optim(par=rep(1, ncol(pure)), fn=sumLossSL, gr=NULL, i, vstSL, shift, power,
            method="L-BFGS-B", lower=0, upper=100)$par
    })
  }

  mix <- do.call(rbind, res)
  mix <- mix / rowSums(mix)
  colnames(mix) <- colnames(pure)
  rownames(mix) <- colnames(x)

  if (format == "matrix") {
    return(mix)
  } else {
    fitted <- pure %*% t(mix)
    cor <- if (missing(shift)) {
      cor(vst(x, alpha), vst(fitted, alpha))
    } else {
      cor(vstSL(x, shift), vstSL(fitted, shift))
    }
    return(list(mix=mix, cor=diag(cor), fitted=fitted))
  }
  
}

#' Collapse technical replicates in a RangedSummarizedExperiment or DESeqDataSet
#'
#' Collapses the columns in \code{object} by summing within levels
#' of a grouping factor \code{groupby}. The purpose of this function
#' is to sum up read counts from technical replicates to create an object
#' with a single column of read counts for each sample.
#' Note: by "technical replicates", we mean multiple sequencing runs of the same
#' library, in constrast to "biological replicates" in which multiple
#' libraries are prepared from separate biological units.
#' Optionally renames the columns of returned object with the levels of the
#' grouping factor.
#' Note: this function is written very simply and
#' can be easily altered to produce other behavior by examining the source code.
#'
#' @param object A \code{RangedSummarizedExperiment} or \code{DESeqDataSet}
#' @param groupby a grouping factor, as long as the columns of object
#' @param run optional, the names of each unique column in object. if provided,
#' a new column \code{runsCollapsed} will be added to the \code{colData}
#' which pastes together the names of \code{run}
#' @param renameCols whether to rename the columns of the returned object
#' using the levels of the grouping factor
#'
#' @return the \code{object} with as many columns as levels in \code{groupby}.
#' This object has assay/count data which is summed from the various
#' columns which are grouped together, and the \code{colData} is subset using
#' the first column for each group in \code{groupby}.
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(m=12)
#' 
#' # make data with two technical replicates for three samples
#' dds$sample <- factor(sample(paste0("sample",rep(1:9, c(2,1,1,2,1,1,2,1,1)))))
#' dds$run <- paste0("run",1:12)
#'
#' ddsColl <- collapseReplicates(dds, dds$sample, dds$run)
#'
#' # examine the colData and column names of the collapsed data
#' colData(ddsColl)
#' colnames(ddsColl)
#'
#' # check that the sum of the counts for "sample1" is the same
#' # as the counts in the "sample1" column in ddsColl
#' matchFirstLevel <- dds$sample == levels(dds$sample)[1]
#' stopifnot(all(rowSums(counts(dds[,matchFirstLevel])) == counts(ddsColl[,1])))
#' 
#' @export
collapseReplicates <- function(object, groupby, run, renameCols=TRUE) {
  if (!is.factor(groupby)) groupby <- factor(groupby)
  groupby <- droplevels(groupby)
  stopifnot(length(groupby) == ncol(object))
  sp <- split(seq(along=groupby), groupby)
  countdata <- sapply(sp, function(i) rowSums(assay(object)[,i,drop=FALSE]))
  mode(countdata) <- "integer"
  colsToKeep <- sapply(sp, `[`, 1)
  collapsed <- object[,colsToKeep]
  dimnames(countdata) <- dimnames(collapsed)
  assay(collapsed) <- countdata
  if (!missing(run)) {
    stopifnot(length(groupby) == length(run))
    colData(collapsed)$runsCollapsed <- sapply(sp, function(i) paste(run[i],collapse=","))
  }
  if (renameCols) {
    colnames(collapsed) <- levels(groupby)
  }
  stopifnot(sum(as.numeric(assay(object))) == sum(as.numeric(assay(collapsed))))
  collapsed
}

#' FPKM: fragments per kilobase per million mapped fragments
#'
#' The following function returns fragment counts normalized
#' per kilobase of feature length per million mapped fragments
#' (by default using a robust estimate of the library size,
#' as in \code{\link{estimateSizeFactors}}).
#'
#' The length of the features (e.g. genes) is calculated one of two ways:
#' (1) If there is a matrix named "avgTxLength" in \code{assays(dds)},
#' this will take precedence in the length normalization.
#' This occurs when using the tximport-DESeq2 pipeline.
#' (2) Otherwise, feature length is calculated 
#' from the \code{rowRanges} of the dds object,
#' if a column \code{basepairs} is not present in \code{mcols(dds)}.
#' The calculated length is the number of basepairs in the union of all \code{GRanges}
#' assigned to a given row of \code{object}, e.g., 
#' the union of all basepairs of exons of a given gene.
#' Note that the second approach over-estimates the gene length
#' (average transcript length, weighted by abundance is a more appropriate
#' normalization for gene counts), and so the FPKM will be an underestimate of the true value.
#' 
#' Note that, when the read/fragment counting has inter-feature dependencies, a strict
#' normalization would not incorporate the basepairs of a feature which
#' overlap another feature. This inter-feature dependence is not taken into
#' consideration in the internal union basepair calculation.
#'
#' @param object a \code{DESeqDataSet}
#' @param robust whether to use size factors to normalize
#' rather than taking the column sums of the raw counts,
#' using the \code{\link{fpm}} function.
#'
#' @return a matrix which is normalized per kilobase of the
#' union of basepairs in the \code{GRangesList} or \code{GRanges}
#' of the mcols(object), and per million of mapped fragments,
#' either using the robust median ratio method (robust=TRUE, default)
#' or using raw counts (robust=FALSE).
#' Defining a column \code{mcols(object)$basepairs} takes
#' precedence over internal calculation of the kilobases for each row.
#'
#' @examples
#'
#' # create a matrix with 1 million counts for the
#' # 2nd and 3rd column, the 1st and 4th have
#' # half and double the counts, respectively.
#' m <- matrix(1e6 * rep(c(.125, .25, .25, .5), each=4),
#'             ncol=4, dimnames=list(1:4,1:4))
#' mode(m) <- "integer"
#' se <- SummarizedExperiment(list(counts=m), colData=DataFrame(sample=1:4))
#' dds <- DESeqDataSet(se, ~ 1)
#' 
#' # create 4 GRanges with lengths: 1, 1, 2, 2.5 Kb
#' gr1 <- GRanges("chr1",IRanges(1,1000)) # 1kb
#' gr2 <- GRanges("chr1",IRanges(c(1,1001),c( 500,1500))) # 1kb
#' gr3 <- GRanges("chr1",IRanges(c(1,1001),c(1000,2000))) # 2kb
#' gr4 <- GRanges("chr1",IRanges(c(1,1001),c(200,1300))) # 500bp
#' rowRanges(dds) <- GRangesList(gr1,gr2,gr3,gr4)
#' 
#' # the raw counts
#' counts(dds)
#'
#' # the FPM values
#' fpm(dds)
#' 
#' # the FPKM values
#' fpkm(dds)
#' 
#' @seealso \code{\link{fpm}}
#'
#' @docType methods
#' @name fpkm
#' @rdname fpkm
#' 
#' @export
fpkm <- function(object, robust=TRUE) {
  fpm <- fpm(object, robust=robust)
  if ("avgTxLength" %in% assayNames(object)) {
    exprs <- 1e3 * fpm / assays(object)[["avgTxLength"]]
    if (robust) {
      sf <- estimateSizeFactorsForMatrix(exprs)
      exprs <- t(t(exprs)/sf)
      return(exprs)
    } else {
      return(exprs)
    }
  }
  if (is.null(mcols(object)$basepairs)) {
    if (is(rowRanges(object), "GRangesList")) {
      ubp <- DataFrame(basepairs = sum(width(reduce(rowRanges(object)))))
    } else if (is(rowRanges(object), "GRanges")) {
      ubp <- DataFrame(basepairs = width(rowRanges(object)))
    }
    if (all(ubp$basepairs == 0)) {
      stop("rowRanges(object) has all ranges of zero width.
the user should instead supply a column, mcols(object)$basepairs,
which will be used to produce FPKM values")
    }
    if (is.null(mcols(mcols(object)))) {
      mcols(object) <- ubp
    } else {
      mcols(ubp) <- DataFrame(type="intermediate",
                              description="count of basepairs in the union of all ranges")
      mcols(object) <- cbind(mcols(object), ubp)
    }
  }
  1e3 * sweep(fpm, 1, mcols(object)$basepairs, "/")
}

#' FPM: fragments per million mapped fragments
#'
#' Calculates either a robust version (default)
#' or the traditional matrix of fragments/counts per million mapped
#' fragments (FPM/CPM).
#' Note: this function is written very simply and
#' can be easily altered to produce other behavior by examining the source code.
#' 
#' @param object a \code{DESeqDataSet}
#' @param robust whether to use size factors to normalize
#' rather than taking the column sums of the raw counts.
#' If TRUE, the size factors and the geometric mean of
#' column sums are multiplied to create a robust library size estimate.
#' Robust normalization is not used if average transcript lengths are present.
#' 
#' @return a matrix which is normalized per million of mapped fragments,
#' either using the robust median ratio method (robust=TRUE, default)
#' or using raw counts (robust=FALSE).
#'
#' @examples
#'
#' # generate a dataset with size factors: .5, 1, 1, 2
#' dds <- makeExampleDESeqDataSet(m = 4, n = 1000,
#'                                interceptMean=log2(1e3),
#'                                interceptSD=0,
#'                                sizeFactors=c(.5,1,1,2),
#'                                dispMeanRel=function(x) .01)
#'
#' # add a few rows with very high count
#' counts(dds)[4:10,] <- 2e5L
#'
#' # in this robust version, the counts are comparable across samples
#' round(head(fpm(dds), 3))
#'
#' # in this column sum version, the counts are still skewed:
#' # sample1 < sample2 & 3 < sample 4
#' round(head(fpm(dds, robust=FALSE), 3))
#'
#' # the column sums of the robust version
#' # are not equal to 1e6, but the
#' # column sums of the non-robust version
#' # are equal to 1e6 by definition
#' 
#' colSums(fpm(dds))/1e6
#' colSums(fpm(dds, robust=FALSE))/1e6
#'
#' @seealso \code{\link{fpkm}}
#'
#' @docType methods
#' @name fpm
#' @rdname fpm
#' 
#' @export
fpm <- function(object, robust=TRUE) {
  # we do something different if average tx lengths are present
  noAvgTxLen <- !("avgTxLength" %in% assayNames(object))
  if (robust & is.null(sizeFactors(object)) & noAvgTxLen) {
    object <- estimateSizeFactors(object)
  }
  k <- counts(object)
  library.sizes <- if (robust & noAvgTxLen) {
    sizeFactors(object) * exp(mean(log(colSums(k))))
  } else {
    colSums(k)
  }
  1e6 * sweep(k, 2, library.sizes, "/")  
}



#' Normalize for gene length
#'
#' Normalize for gene length using the output of transcript abundance estimators
#'
#' This function is deprecated and moved to a new general purpose package,
#' tximport, which will be added to Bioconductor.
#'
#' @param ... ...
#' 
#' @export
normalizeGeneLength <- function(...) {
  .Deprecated("tximport, a separate package on Bioconductor")
}

#' Normalized counts transformation
#'
#' A simple function for creating a \code{\link{DESeqTransform}}
#' object after applying: \code{f(count(dds,normalized=TRUE) + pc)}.
#' 
#' @param object a DESeqDataSet object
#' @param f a function to apply to normalized counts
#' @param pc a pseudocount to add to normalized counts
#' 
#' @seealso \code{\link{varianceStabilizingTransformation}}, \code{\link{rlog}}
#' 
#' @export
normTransform <- function(object, f=log2, pc=1) {
  if (is.null(colnames(object))) {
    colnames(object) <- seq_len(ncol(object))
  }
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  nt <- f(counts(object, normalized=TRUE) + pc)
  se <- SummarizedExperiment(
    assays = nt,
    colData = colData(object),
    rowRanges = rowRanges(object),
    metadata = metadata(object))
  DESeqTransform(se)
}
