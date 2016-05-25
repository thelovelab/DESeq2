#' Collapse technical replicates in a RangedSummarizedExperiment or DESeqDataSet
#'
#' Collapses the columns in \code{object} by summing within levels
#' of a grouping factor \code{groupby}. The purpose of this function
#' is to sum up read counts from technical replicates to create an object
#' with a single column of read counts for each sample.
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
    if (class(rowRanges(object)) == "GRangesList") {
      ubp <- DataFrame(basepairs = sum(width(reduce(rowRanges(object)))))
    } else if (class(rowRanges(object)) == "GRanges") {
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
#' object after applying: f(count + pc).
#' 
#' @param object a DESeqDataSet object
#' @param f a function to apply to normalized counts
#' @param pc a pseudocount to add to normalized counts
#' 
#' @seealso \code{\link{varianceStabilizingTransformation}}, \code{\link{rlog}}
#' 
#' @export
normTransform <- function(object, f=log2, pc=1) {
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


#####################
# unexported
#####################









# function to split up DESeqDataSet by rows during easily parallelizable steps

# TODO: recombining the resulting DESeqDataSets using rbind() is a bit wasteful,
# as the count matrix and GRanges from the original object are unchanged.

DESeqParallel <- function(object, test, fitType, betaPrior, full, reduced, quiet, modelMatrix, modelMatrixType, BPPARAM) {

  # size factors already estimated or supplied
  # break up the object into equal sized chunks
  # to be fed to the different workers
  object <- getBaseMeansAndVariances(object)
  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]
  nworkers <- BPPARAM$workers
  idx <- factor(sort(rep(seq_len(nworkers),length=nrow(objectNZ))))

  if (missing(modelMatrixType)) {
    modelMatrixType <- NULL
  }

  # if no reps, treat samples as replicates and print warning
  noReps <- checkForExperimentalReplicates(object, modelMatrix)
  if (noReps) {
    designIn <- design(objectNZ)
    design(objectNZ) <- formula(~ 1)
  }
  
  # first parallel execution: gene-wise dispersion estimates
  if (!quiet) message("estimating dispersions")
  if (!quiet) message(paste("gene-wise dispersion estimates:",nworkers,"workers"))
  objectNZ <- do.call(rbind, bplapply(levels(idx), function(l) {
    estimateDispersionsGeneEst(objectNZ[idx == l,,drop=FALSE], quiet=TRUE, modelMatrix=modelMatrix)
  }, BPPARAM=BPPARAM))

  # the dispersion fit and dispersion prior are estimated over all rows
  if (!quiet) message("mean-dispersion relationship") 
  objectNZ <- estimateDispersionsFit(objectNZ, fitType=fitType)
  dispPriorVar <- estimateDispersionsPriorVar(objectNZ, modelMatrix=modelMatrix)

  # need to condition on whether a beta prior needs to be fit
  if (betaPrior) {

    # if so,

    # second parallel execution: fit the final dispersion estimates and MLE betas 
    if (!quiet) message(paste("final dispersion estimates, MLE betas:",nworkers,"workers"))
    objectNZ <- do.call(rbind, bplapply(levels(idx), function(l) {
      objectNZSub <- estimateDispersionsMAP(objectNZ[idx == l,,drop=FALSE],
                                            dispPriorVar=dispPriorVar, quiet=TRUE)
      # replace design
      if (noReps) design(objectNZSub) <- designIn
      estimateMLEForBetaPriorVar(objectNZSub)
    }, BPPARAM=BPPARAM))

    # replace design
    if (noReps) design(objectNZ) <- designIn
    
    # need to set standard model matrix for LRT with beta prior
    if (test == "LRT") {
      attr(object, "modelMatrixType") <- "standard"
      attr(objectNZ, "modelMatrixType") <- "standard"
      modelMatrixType <- "standard"
    }
    
    # the beta prior is estimated over all rows
    betaPriorVar <- estimateBetaPriorVar(objectNZ)
    
    # the third parallel execution: the final GLM and statistics
    if (!quiet) message(paste("fitting model and testing:",nworkers,"workers"))
    if (test == "Wald") {

      objectNZ <- do.call(rbind, bplapply(levels(idx), function(l) {
        nbinomWaldTest(objectNZ[idx == l,,drop=FALSE], betaPriorVar=betaPriorVar,
                       quiet=TRUE, modelMatrixType=modelMatrixType)
      }, BPPARAM=BPPARAM))
    } else if (test == "LRT") {
      objectNZ <- do.call(rbind, bplapply(levels(idx), function(l) {
        nbinomLRT(objectNZ[idx == l,,drop=FALSE], full=full, reduced=reduced,
                  betaPrior=betaPrior, betaPriorVar=betaPriorVar, quiet=TRUE)
      }, BPPARAM=BPPARAM))
    }
    
  } else {
    
    # or, if no beta prior to fit,
 
    # second parallel execution: fit the final dispersion estimates and the final GLM and statistics
    if (!quiet) message(paste("final dispersion estimates, fitting model and testing:",nworkers,"workers"))
    if (test == "Wald") {
      objectNZ <- do.call(rbind, bplapply(levels(idx), function(l) {
        objectNZSub <- estimateDispersionsMAP(objectNZ[idx == l,,drop=FALSE],
                                              dispPriorVar=dispPriorVar, quiet=TRUE, modelMatrix=modelMatrix)
        # replace design
        if (noReps) design(objectNZSub) <- designIn
        nbinomWaldTest(objectNZSub, betaPrior=betaPrior,
                       quiet=TRUE, modelMatrix=modelMatrix, modelMatrixType="standard")
      }, BPPARAM=BPPARAM))
    } else if (test == "LRT") {
      objectNZ <- do.call(rbind, bplapply(levels(idx), function(l) {
        objectNZSub <- estimateDispersionsMAP(objectNZ[idx == l,,drop=FALSE],
                                              dispPriorVar=dispPriorVar, quiet=TRUE, modelMatrix=modelMatrix)
        # replace design
        if (noReps) design(objectNZSub) <- designIn
        nbinomLRT(objectNZSub, full=full, reduced=reduced, quiet=TRUE)
      }, BPPARAM=BPPARAM))
    } 
  }

  outMcols <- buildDataFrameWithNARows(mcols(objectNZ), mcols(object)$allZero)
  mcols(outMcols) <- mcols(mcols(objectNZ))
  outMu <- buildMatrixWithNARows(assays(objectNZ)[["mu"]], mcols(object)$allZero)
  outCooks <- buildMatrixWithNARows(assays(objectNZ)[["cooks"]], mcols(object)$allZero)

  # now backfill any columns in rowRanges which existed before running DESeq()
  # and which are not of type "intermediate" or "results"
  object <- sanitizeRowRanges(object)
  inMcols <- mcols(object)
  namesCols <- names(mcols(object))
  inputCols <- namesCols[! mcols(mcols(object))$type %in% c("intermediate","results")]
  for (var in inputCols) {
    outMcols[var] <- inMcols[var]
  }
  
  mcols(object) <- outMcols
  object <- getBaseMeansAndVariances(object)
  assays(object)[["mu"]] <- outMu
  assays(object)[["cooks"]] <- outCooks

  attrNames <- c("dispersionFunction","modelMatrixType","betaPrior",
                 "betaPriorVar","modelMatrix","test","dispModelMatrix")
  for (attrName in attrNames) {
    attr(object, attrName) <- attr(objectNZ, attrName)
  }
  
  object
}
