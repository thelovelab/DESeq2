context("outlier")
test_that("outlier filtering and replacement works as expected", {
  # test filtering and replacement
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100, m=12, dispMeanRel = function(x) 4/x + .5)
  counts(dds)[1,] <- rep(0L, 12)
  counts(dds)[2,] <- c(100000L, rep(10L, 11))
  counts(dds)[3,] <- c(100000L, rep(0L, 11))
  dds0 <- DESeq(dds, minReplicatesForReplace=Inf)
  dds1 <- DESeq(dds, minReplicatesForReplace=6)
  pval0 <- results(dds0)[1:3,"pvalue"]
  pval <- results(dds1)[1:3,"pvalue"]
  LFC0 <- results(dds0)[1:3,"log2FoldChange"]
  LFC <- results(dds1)[1:3,"log2FoldChange"]

  # filtered
  expect_true(all(is.na(pval0)))
  # not filtered
  expect_true(all(!is.na(pval[2:3])))
  # counts still the same
  expect_true(all(counts(dds1)==counts(dds)))
  # first is NA
  expect_true(is.na(LFC[1]))
  # replaced, reduced LFC
  expect_true(abs(LFC[2]) < abs(LFC0[2]))
  # replaced, LFC now zero
  expect_true(LFC[3] == 0)
  idx <- which(!mcols(dds1)$replace)
  # the pvalue for those not replaced is equal
  expect_equal(results(dds1)$pvalue[idx], results(dds0)$pvalue[idx])

  # check that outlier filtering catches throughout range of mu
  beta0 <- seq(from=1,to=16,length=100)
  idx <- rep(rep(c(TRUE,FALSE),c(1,9)),10)
  set.seed(1)
  #par(mfrow=c(2,3))
  for (disp0 in c(.01,.1)) {
    for (m in c(10,20,80)) {
      dds <- makeExampleDESeqDataSet(n=100, m=m, interceptMean=beta0, interceptSD=0,
                                     dispMeanRel=function(x) disp0)
      counts(dds)[idx,1] <- as.integer(1000 * 2^beta0[idx])
      dds <- DESeq(dds, minReplicatesForReplace=Inf, quiet=TRUE, fitType="mean")
      res <- results(dds)
      cutoff <- qf(.99, 2, m-2)
      outlierCooks <- assays(dds)[["cooks"]][idx,1] > cutoff
      nonoutlierCooks <- mcols(dds)$maxCooks[!idx] < cutoff
      expect_true(all(is.na(res$pvalue[idx])))
      expect_true(all(outlierCooks))
      expect_true(all(nonoutlierCooks))
      col <- rep("black", 100)
      col[idx] <- "blue"
      #plot(2^beta0, mcols(dds)$maxCooks, col=col, log="xy",
      #     main=paste(m,"-",disp0), xlab="mean", ylab="cooks")
      #abline(h=qf(.99,2,m-2))
    }
  }

  dds <- makeExampleDESeqDataSet(n=100)
  counts(dds)[1,1] <- 1000000L
  dds <- DESeq(dds, test="LRT", reduced=~1, minReplicatesForReplace=6)

  # test replace function
  dds <- makeExampleDESeqDataSet(n=100,m=4)
  expect_error(replaceOutliers(dds))
  dds <- DESeq(dds)
  expect_error(replaceOutliers(dds, minReplicates=2))

  # check model matrix standard bug
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100, m=20)
  counts(dds)[1,] <- c(100000L, rep(0L, 19))
  dds <- DESeq(dds, modelMatrixType="standard")
})

test_that("outlier filtering doesn't flag small counts", {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100, m=8, dispMeanRel=function(x) 0.01)
  counts(dds)[1,] <- c(0L, 0L, 0L, 100L, 2100L, 2200L, 2300L, 2400L)
  counts(dds)[2:3,1] <- 100000L
  counts(dds)[4,] <- rep(0L, 8)
  dds <- DESeq(dds, fitType="mean")
  res <- results(dds)
  expect_true(!is.na(res$pvalue[1]))
  expect_true(all(is.na(res$pvalue[2:3])))
})
