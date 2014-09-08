test_parallel <- function() {
  dispMeanRel <- function(x) (4/x + .1) * exp(rnorm(1000,0,sqrt(.5)))
  set.seed(1)
  dds0 <- makeExampleDESeqDataSet(dispMeanRel=dispMeanRel)

  # the following is an example of a simple parallelizable DESeq() run
  # without outlier replacement. see DESeq2:::DESeqParallel for the code
  # which is actually used in DESeq()
  
  parallel <- 6
  idx <- factor(sort(rep(seq_len(parallel),length=nrow(dds0))))

  ### BEGINNING ###
  
  dds <- estimateSizeFactors(dds0)
  dds <- do.call(rbind, lapply(levels(idx), function(l) {
    estimateDispersionsGeneEst(dds[idx == l,,drop=FALSE])
  }))
  dds <- estimateDispersionsFit(dds)
  dispPriorVar <- estimateDispersionsPriorVar(dds)
  dds <- do.call(rbind, lapply(levels(idx), function(l) {
    ddsSub <- estimateDispersionsMAP(dds[idx == l,,drop=FALSE], dispPriorVar=dispPriorVar)
    estimateMLEForBetaPriorVar(ddsSub)
  }))
  betaPriorVar <- estimateBetaPriorVar(dds)
  dds <- do.call(rbind, lapply(levels(idx), function(l) {
    nbinomWaldTest(dds[idx == l,,drop=FALSE], betaPriorVar=betaPriorVar)
  }))  

  ### END ###
  
  res1 <- results(dds)
  
  dds2 <- DESeq(dds0)
  res2 <- results(dds2)

  checkEqualsNumeric(mcols(dds)$dispGeneEst, mcols(dds2)$dispGeneEst)
  checkEqualsNumeric(mcols(dds)$dispFit, mcols(dds2)$dispFit)
  checkEqualsNumeric(mcols(dds)$dispMAP, mcols(dds2)$dispMAP)
  checkEqualsNumeric(mcols(dds)$dispersion, mcols(dds2)$dispersion)
  checkEqualsNumeric(attr(dispersionFunction(dds), "dispPriorVar"),
                     attr(dispersionFunction(dds2), "dispPriorVar"))
  checkEqualsNumeric(attr(dispersionFunction(dds), "varLogDispEsts"),
                     attr(dispersionFunction(dds2), "varLogDispEsts"))
  checkEqualsNumeric(mcols(dds)$MLE_condition_B_vs_A,
                     mcols(dds2)$MLE_condition_B_vs_A)
  checkEqualsNumeric(attr(dds, "betaPriorVar"),
                     attr(dds2, "betaPriorVar"))
  checkEqualsNumeric(mcols(dds)$conditionB, mcols(dds2)$conditionB)  
  checkEqualsNumeric(res1$log2FoldChange, res2$log2FoldChange)
  checkEqualsNumeric(res1$pvalue, res2$pvalue)

  library(BiocParallel)
  register(SerialParam())
  dds <- DESeq(dds0, parallel=6)
  
}
