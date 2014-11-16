#' Collapse technical replicates in a SummarizedExperiment or DESeqDataSet
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
#' @param object A \code{SummarizedExperiment} or \code{DESeqDataSet}
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
  assay(collapsed) <- countdata
  if (!missing(run)) {
    stopifnot(length(groupby) == length(run))
    colData(collapsed)$runsCollapsed <- sapply(sp, function(i) paste(run[i],collapse=","))
  }
  if (renameCols) {
    colnames(collapsed) <- levels(groupby)
  }
  stopifnot(sum(assay(object)) == sum(assay(collapsed)))
  collapsed
}

#' FPKM: fragments per kilobase per million mapped fragments
#'
#' The following function returns fragment counts normalized
#' per kilobase of feature length per million mapped fragments
#' (by default using a robust estimate of the library size,
#' as in \code{\link{estimateSizeFactors}}).
#'
#' Note: the kilobase length of the features is calculated from the \code{rowData}
#' if a column \code{basepairs} is not present in \code{mcols(dds)}.
#' This is the number of basepairs in the union of all \code{GRanges}
#' assigned to a given row of \code{object}, e.g., 
#' the union of all basepairs of exons of a given gene.
#' When the read/fragment counting is interfeature dependent, a strict
#' normalization would not incorporate the basepairs of a feature which
#' overlap another feature. This interfeature dependence is not taken into
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
#' se <- SummarizedExperiment(m, colData=DataFrame(sample=1:4))
#' dds <- DESeqDataSet(se, ~ 1)
#' 
#' # create 4 GRanges with lengths: 1, 1, 2, 2.5 Kb
#' gr1 <- GRanges("chr1",IRanges(1,1000))
#' gr2 <- GRanges("chr1",IRanges(c(1,501),c(500,1000)))
#' gr3 <- GRanges("chr1",IRanges(c(1,1001),c(1000,2000)))
#' gr4 <- GRanges("chr1",IRanges(c(1,1001,2001),c(500,3000,3000)))
#' rowData(dds) <- GRangesList(gr1,gr2,gr3,gr4)
#' 
#' # the raw counts
#' counts(dds)
#' 
#' # the FPKM values
#' fpkm(dds)
#' 
#' # held constant per 1 million fragments
#' counts(dds) <- counts(dds) * 2L
#' round(fpkm(dds))
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
  if (is.null(mcols(object)$basepairs)) {
    if (class(rowData(object)) == "GRangesList") {
      ubp <- DataFrame(basepairs = sum(width(reduce(rowData(object)))))
    } else if (class(rowData(object)) == "GRanges") {
      ubp <- DataFrame(basepairs = width(rowData(object)))
    }
    if (all(ubp$basepairs == 0)) {
      stop("rowData(object) has all ranges of zero width.
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
#' 
#' @return a matrix which is normalized per million of mapped fragments,
#' either using the robust median ratio method (robust=TRUE, default)
#' or using raw counts (robust=FALSE).
#'
#' @examples
#'
#' # generate a dataset with size factors: .5, 1, 1, 2
#' dds <- makeExampleDESeqDataSet(m = 4, n = 1e3,
#'                                interceptMean=log2(1e3),
#'                                interceptSD=0,
#'                                sizeFactors=c(.5,1,1,2),
#'                                dispMeanRel=function(x) .01)
#'
#' # examine the column sums
#' # then add 1 million counts over rows 11:15 for each sample
#' colSums(counts(dds))/1e6
#' counts(dds)[11:15,] <- 2e5L
#' colSums(counts(dds))/1e6
#'
#' # the robust FPM treats the samples
#' # relatively equally
#' head(fpm(dds), 3)
#'
#' # the non-robust version is thrown
#' # off by the 5 rows with large counts
#' head(fpm(dds, robust=FALSE), 3)
#'
#' # the column sums of the robust version
#' # are not equal to 1e6, but the
#' # column sums of the non-robust version
#' # are equal to 1e6 by definition
#' colSums(fpm(dds))/1e6
#' colSums(fpm(dds, robust=FALSE))/1e6
#'
#' # the total sum is equal for both methods
#' sum(fpm(dds))
#' sum(fpm(dds, robust=FALSE))
#'
#' @seealso \code{\link{fpkm}}
#'
#' @docType methods
#' @name fpm
#' @rdname fpm
#' 
#' @export
fpm <- function(object, robust=TRUE) {
  if (robust & is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  if (robust) {
    k <- counts(object, normalized=TRUE)
    1e6 * ncol(k) * k / sum(k)
  } else {
    k <- counts(object)
    1e6 * sweep(k, 2, colSums(k), "/")
  }
}



#####################
# unexported
#####################

# function to split up DESeqDataSet by rows during easily parallelizable steps

# TODO: recombining the resulting DESeqDataSets using rbind() is a bit wasteful,
# as the count matrix and GRanges from the original object are unchanged.

DESeqParallel <- function(object, test, fitType, betaPrior, full, reduced, quiet, modelMatrixType, BPPARAM) {

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
  
  # first parallel execution: gene-wise dispersion estimates
  if (!quiet) message("estimating dispersions")
  if (!quiet) message(paste("gene-wise dispersion estimates:",nworkers,"workers"))
  objectNZ <- do.call(rbind, bplapply(levels(idx), function(l) {
    estimateDispersionsGeneEst(objectNZ[idx == l,,drop=FALSE], quiet=TRUE)
  }, BPPARAM=BPPARAM))

  # the dispersion fit and dispersion prior are estimated over all rows
  if (!quiet) message("mean-dispersion relationship") 
  objectNZ <- estimateDispersionsFit(objectNZ, fitType=fitType)
  dispPriorVar <- estimateDispersionsPriorVar(objectNZ)

  # need to condition on whether a beta prior needs to be fit
  if (betaPrior) {

    # if so,

    # second parallel execution: fit the final dispersion estimates and MLE betas 
    if (!quiet) message(paste("final dispersion estimates, MLE betas:",nworkers,"workers"))
    objectNZ <- do.call(rbind, bplapply(levels(idx), function(l) {
      objectNZSub <- estimateDispersionsMAP(objectNZ[idx == l,,drop=FALSE],
                                            dispPriorVar=dispPriorVar, quiet=TRUE)
      estimateMLEForBetaPriorVar(objectNZSub)
    }, BPPARAM=BPPARAM))

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
      if (missing(modelMatrixType)) {
        modelMatrixType <- "standard"
      }
      objectNZ <- do.call(rbind, bplapply(levels(idx), function(l) {
        nbinomLRT(objectNZ[idx == l,,drop=FALSE], full=full, reduced=reduced,
                  betaPrior=betaPrior, betaPriorVar=betaPriorVar,
                  quiet=TRUE, modelMatrixType=modelMatrixType)
      }, BPPARAM=BPPARAM))
    }
    
  } else {
    
    # or, if no beta prior to fit,
 
    # second parallel execution: fit the final dispersion estimates and the final GLM and statistics
    if (!quiet) message(paste("final dispersion estimates, fitting model and testing:",nworkers,"workers"))
    if (test == "Wald") {
      objectNZ <- do.call(rbind, bplapply(levels(idx), function(l) {
        objectNZSub <- estimateDispersionsMAP(objectNZ[idx == l,,drop=FALSE],
                                              dispPriorVar=dispPriorVar, quiet=TRUE)
        nbinomWaldTest(objectNZSub, betaPrior=betaPrior,
                       quiet=TRUE, modelMatrixType="standard")
      }, BPPARAM=BPPARAM))
    } else if (test == "LRT") {
      objectNZ <- do.call(rbind, bplapply(levels(idx), function(l) {
        objectNZSub <- estimateDispersionsMAP(objectNZ[idx == l,,drop=FALSE],
                                              dispPriorVar=dispPriorVar, quiet=TRUE)
        nbinomLRT(objectNZSub, full=full, reduced=reduced,
                  quiet=TRUE, modelMatrixType="standard")
      }, BPPARAM=BPPARAM))
    }
    
  }

  outMcols <- buildDataFrameWithNARows(mcols(objectNZ), mcols(object)$allZero)
  mcols(outMcols) <- mcols(mcols(objectNZ))
  outMu <- buildMatrixWithNARows(assays(objectNZ)[["mu"]], mcols(object)$allZero)
  outCooks <- buildMatrixWithNARows(assays(objectNZ)[["cooks"]], mcols(object)$allZero)
  
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
