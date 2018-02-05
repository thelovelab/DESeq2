context("nbinomWald")
test_that("nbinomWald throws various errors and works with edge cases",{
  dds <- makeExampleDESeqDataSet(n=100, m=4)
  expect_error(nbinomWaldTest(dds))
  expect_error(nbinomLRT(dds))
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  mm <- model.matrix(~ condition, colData(dds))
  mm0 <- model.matrix(~ 1, colData(dds))
  expect_error(nbinomWaldTest(dds, betaPrior=TRUE, modelMatrix=mm))
  expect_error(nbinomLRT(dds, betaPrior=TRUE, full=mm, reduced=mm0))
  expect_error(nbinomWaldTest(dds, betaPrior=FALSE, modelMatrixType="expanded"))
  expect_error(nbinomLRT(dds, betaPrior=FALSE, modelMatrixType="expanded"))
  dds2 <- estimateMLEForBetaPriorVar(dds)
  estimateBetaPriorVar(dds2, betaPriorMethod="quantile")
  dds <- nbinomWaldTest(dds, modelMatrixType="standard")
  covarianceMatrix(dds, 1)

  # changing 'df'
  dds <- makeExampleDESeqDataSet(n=100, m=4)
  counts(dds)[1:4,] <- rep(0L, 16)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  round(head(results(dds)$pvalue,8),3)
  dds <- nbinomWaldTest(dds, useT=TRUE, df=rep(1,100))
  round(head(results(dds)$pvalue,8),3)
  
  # try nbinom after no fitted dispersions
  dds <- makeExampleDESeqDataSet(n=100, m=4)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
})

test_that("useT uses proper degrees of freedom", {

  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=200, m=15)
  counts(dds)[101:105,] <- 0L
  dds$condition <- factor(rep(c("A","B","C"),each=5))
  dds <- DESeq(dds, useT=TRUE)
  dds <- removeResults(dds)
  w <- matrix(1, nrow=nrow(dds), ncol=ncol(dds))
  w[1:100,1] <- 0
  w[1,c(1:4,6:9,11:14)] <- 0
  assays(dds)[["weights"]] <- w
  dds <- DESeq(dds, useT=TRUE)
  res <- results(dds)
  expect_true(is.na(res$pvalue[1]))
  expect_true(mcols(dds)$tDegreesFreedom[2] == 15-1-3)
  expect_true(res$pvalue[2] == 2*pt(abs(res$stat[2]), df=15-1-3, lower.tail=FALSE))

  # also lfcThreshold
  res <- results(dds, lfcThreshold=1, altHypothesis="greaterAbs")
  idx <- which(res$log2FoldChange > 1 & !is.na(res$pvalue))[1]
  expect_true(res$pvalue[idx] == 2 * pt(res$stat[idx], df=15-1-3, lower.tail=FALSE))
  #
  res <- results(dds, lfcThreshold=1, altHypothesis="greater")
  idx <- which(res$log2FoldChange > 1 & !is.na(res$pvalue))[1]
  expect_true(res$pvalue[idx] == pt(res$stat[idx], df=15-1-3, lower.tail=FALSE))
  #
  res <- results(dds, lfcThreshold=1, altHypothesis="less")
  idx <- which(res$log2FoldChange < -1 & !is.na(res$pvalue))[1]
  expect_true(res$pvalue[idx] == pt(-1 * res$stat[idx], df=15-1-3, lower.tail=FALSE))
  #
  res <- results(dds, lfcThreshold=1, altHypothesis="lessAbs")
  idx <- which(abs(res$log2FoldChange) < 1 & !is.na(res$pvalue))[1]
  expect_true(res$pvalue[idx] == pt(res$stat[idx], df=15-1-3, lower.tail=FALSE))    
  
  # also novel contrasts
  res <- results(dds, contrast=c("condition","C","B"))
  expect_true(is.na(res$pvalue[1]))
  expect_true(mcols(dds)$tDegreesFreedom[2] == 15-1-3)
  expect_true(res$pvalue[2] == 2*pt(abs(res$stat[2]), df=15-1-3, lower.tail=FALSE))
  
})
