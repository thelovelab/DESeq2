context("parallel")
test_that("parallel execution works as expected", {

  set.seed(1)
  dds0 <- makeExampleDESeqDataSet(n=100)
  counts(dds0)[51:60,] <- 0L

  nworkers <- 4
  idx <- factor(sort(rep(seq_len(nworkers),length=nrow(dds0))))

  dds <- estimateSizeFactors(dds0)
  dds <- do.call(rbind, lapply(levels(idx), function(l) {
    estimateDispersionsGeneEst(dds[idx == l,,drop=FALSE])
  }))
  dds <- estimateDispersionsFit(dds)
  dispPriorVar <- estimateDispersionsPriorVar(dds)
  dds <- do.call(rbind, lapply(levels(idx), function(l) {
    ddsSub <- estimateDispersionsMAP(dds[idx == l,,drop=FALSE], dispPriorVar=dispPriorVar)
    nbinomWaldTest(ddsSub)
  }))

  res1 <- results(dds)

  dds2 <- DESeq(dds0)
  res2 <- results(dds2)

  expect_equal(mcols(dds)$dispGeneEst, mcols(dds2)$dispGeneEst)
  expect_equal(mcols(dds)$dispFit, mcols(dds2)$dispFit)
  expect_equal(mcols(dds)$dispMAP, mcols(dds2)$dispMAP)
  expect_equal(mcols(dds)$dispersion, mcols(dds2)$dispersion)
  expect_equal(attr(dispersionFunction(dds), "dispPriorVar"),
               attr(dispersionFunction(dds2), "dispPriorVar"))
  expect_equal(attr(dispersionFunction(dds), "varLogDispEsts"),
               attr(dispersionFunction(dds2), "varLogDispEsts"))
  expect_equal(mcols(dds)$WaldStatistic_condition_B_vs_A,
               mcols(dds2)$WaldStatistic_condition_B_vs_A)
  expect_equal(res1$pvalue, res2$pvalue)

  # try multicore
  if (FALSE) {
    library("BiocParallel")
    register(MulticoreParam(4))

    # examine metadata proliferation
    metadata(dds0)$foo <- "bar"
    dds3 <- DESeq(dds0, parallel=TRUE)
    metadata(dds3)
    
    expect_true(length(metadata(dds3)) == 2)
    
    res3 <- results(dds3, parallel=TRUE)
    expect_equal(res2$pvalue, res3$pvalue)

    # LRT
    dds.lrt <- DESeq(dds0, parallel=TRUE, test="LRT", reduced=~1)
    
    # lfcShrink parallel test
    dds <- DESeq(dds0)
    # normal
    res <- lfcShrink(dds, coef=2)
    res2 <- lfcShrink(dds, coef=2, parallel=TRUE)
    expect_equal(res$log2FoldChange, res2$log2FoldChange)
    # apeglm
    res <- lfcShrink(dds, coef=2, type="apeglm", svalue=TRUE)
    res2 <- lfcShrink(dds, coef=2, type="apeglm", parallel=TRUE, svalue=TRUE)
    expect_equal(res$log2FoldChange, res2$log2FoldChange)
    expect_equal(res$svalue, res2$svalue)
  }
  
})
