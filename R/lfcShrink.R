#' Shrink log2 fold changes
#'
#' This function adds shrunken log2 fold changes (LFC) to a
#' results table which was run without LFC moderation.
#'
#' At the moment, shrinkage cannot be applied to coefficients
#' in a model with interaction terms, but this will hopefully
#' be added as a feature in the devel cycle of version 1.17.
#' 
#' @param dds a DESeqDataSet object, which has been run through
#' \code{\link{DESeq}}, or at the least, \code{\link{estimateDispersions}}
#' @param coef the number of the coefficient (LFC) to shrink,
#' consult \code{resultsNames(dds)} after running \code{DESeq(dds, betaPrior=FALSE)}.
#' only \code{coef} or \code{contrast} can be specified, not both
#' @param contrast see argument description in \code{\link{results}}.
#' only \code{coef} or \code{contrast} can be specified, not both
#' @param res a DESeqResults object (can be missing)
#' @param type at this time, ignored argument, because only one
#' shrinkage estimator, but more to come!
#'
#' @return if \code{res} is not missing, a DESeqResults object with
#' the \code{log2FoldChange} column replaced with a shrunken LFC.
#' If \code{res} is missing, just the shrunken LFC vector.
#'
#' @export
#' 
#' @examples
#'
#'  dds <- makeExampleDESeqDataSet(betaSD=1)
#'  dds <- DESeq(dds, betaPrior=FALSE)
#'  res <- results(dds)
#'  res.shr <- lfcShrink(dds=dds, coef=2, res=res)
#'  res.shr <- lfcShrink(dds=dds, contrast=c("condition","B","A"), res=res)
#' 
lfcShrink <- function(dds, coef, contrast, res, type="normal") {
  if (is.null(dispersions(dds))) {
    stop("lfcShrink requires dispersion estimates, first call estimateDispersions()")
  }

  # TODO: lfcThreshold
  
  # match the shrinkage type
  type <- match.arg(type, choices=c("normal"))

  # fit MLE coefficients... TODO skip this step
  dds <- estimateMLEForBetaPriorVar(dds)

  stopifnot(missing(coef) | missing(contrast))
  if (missing(contrast)) {
    modelMatrixType <- "standard"
  } else {
    modelMatrixType <- "expanded"
  }
  attr(dds,"modelMatrixType") <- modelMatrixType
  betaPriorVar <- estimateBetaPriorVar(dds)

  dds.shr <- nbinomWaldTest(dds,
                            betaPrior=TRUE,
                            betaPriorVar=betaPriorVar,
                            modelMatrixType=modelMatrixType,
                            quiet=TRUE)

  if (missing(contrast)) {
    rn <- resultsNames(dds.shr)
    res.shr <- results(dds.shr, name=rn[coef])
  } else {
    res.shr <- results(dds.shr, contrast=contrast)
  }

  # if a DESeqResults object was provided
  if (!missing(res)) {
    res$log2FoldChange <- res.shr$log2FoldChange
    res$lfcSE <- res.shr$lfcSE
    mcols(res)$description[2:3] <- mcols(res.shr)$description[2:3]
    deseq2.version <- packageVersion("DESeq2")
    priorInfo(res) <- list(type="normal", package="DESeq2", version=deseq2.version,
                           betaPriorVar=betaPriorVar)
    return(res)
  } else {
    return(res.shr$log2FoldChange)
  }
}
