# function to split up DESeqDataSet by rows during easily parallelizable steps

# TODO: recombining the resulting DESeqDataSets using rbind() is a bit wasteful,
# as the count matrix and GRanges from the original object are unchanged

DESeqParallel <- function(object, test, fitType, betaPrior, full, reduced,
                          quiet, modelMatrix, useT, minmu, BPPARAM) {

  nworkers <- BPPARAM$workers
  idx <- factor(sort(rep(seq_len(nworkers),length.out=nrow(object))))

  checkForExperimentalReplicates(object, modelMatrix)
  
  # first parallel execution: gene-wise dispersion estimates
  if (!quiet) message("estimating dispersions")
  if (!quiet) message(paste("gene-wise dispersion estimates:",nworkers,"workers"))

  object <- do.call(rbind, bplapply(levels(idx), function(l) {
    estimateDispersionsGeneEst(object[idx == l,], quiet=TRUE, modelMatrix=modelMatrix, minmu=minmu)
  }, BPPARAM=BPPARAM))

  # the dispersion fit and dispersion prior are estimated over all rows
  if (!quiet) message("mean-dispersion relationship") 
  object <- estimateDispersionsFit(object, fitType=fitType)
  dispPriorVar <- estimateDispersionsPriorVar(object, modelMatrix=modelMatrix)

  # need to condition on whether a beta prior needs to be fit
  if (betaPrior) {
    # second parallel execution: fit the final dispersion estimates and MLE betas 
    if (!quiet) message(paste("final dispersion estimates, MLE betas:",nworkers,"workers"))
    object <- do.call(rbind, bplapply(levels(idx), function(l) {
      objectSub <- estimateDispersionsMAP(object[idx == l,],
                                          dispPriorVar=dispPriorVar, quiet=TRUE)
      estimateMLEForBetaPriorVar(objectSub)
    }, BPPARAM=BPPARAM))
    # the beta prior is estimated over all rows
    betaPriorVar <- estimateBetaPriorVar(object)
    # the third parallel execution: the final GLM and statistics
    if (!quiet) message(paste("fitting model and testing:",nworkers,"workers"))
    object <- do.call(rbind, bplapply(levels(idx), function(l) {
      nbinomWaldTest(object[idx == l,],
                     betaPrior=TRUE,
                     betaPriorVar=betaPriorVar,
                     quiet=TRUE, useT=useT, minmu=minmu)
    }, BPPARAM=BPPARAM))
  } else {
    # or, if no beta prior to fit,
    # second parallel execution: fit the final dispersion estimates and the final GLM and statistics
    if (!quiet) message(paste("final dispersion estimates, fitting model and testing:",nworkers,"workers"))
    if (test == "Wald") {
      object <- do.call(rbind, bplapply(levels(idx), function(l) {
        objectSub <- estimateDispersionsMAP(object[idx == l,],
                                            dispPriorVar=dispPriorVar, quiet=TRUE, modelMatrix=modelMatrix)
        nbinomWaldTest(objectSub, betaPrior=FALSE,
                       quiet=TRUE, modelMatrix=modelMatrix,
                       useT=useT, minmu=minmu)
      }, BPPARAM=BPPARAM))
    } else if (test == "LRT") {
      object <- do.call(rbind, bplapply(levels(idx), function(l) {
        objectSub <- estimateDispersionsMAP(object[idx == l,],
                                            dispPriorVar=dispPriorVar, quiet=TRUE, modelMatrix=modelMatrix)
        nbinomLRT(objectSub, full=full, reduced=reduced, quiet=TRUE, minmu=minmu)
      }, BPPARAM=BPPARAM))
    } 
  }
  object
}
