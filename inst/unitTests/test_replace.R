test_replace <- function() {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100, m=12, dispMeanRel = function(x) 4/x + .5)
  counts(dds)[1,] <- rep(0L, 12)
  counts(dds)[2,] <- c(100000L, rep(10L, 11))
  counts(dds)[3,] <- c(100000L, rep(0L, 11))
  dds0 <- DESeq(dds, minReplicatesForReplace=Inf)
  dds1 <- DESeq(dds, minReplicatesForReplace=6)

  head(results(dds0),3)
  head(results(dds1),3)

  pval0 <- results(dds0)[1:3,"pvalue"]
  pval <- results(dds1)[1:3,"pvalue"]
  LFC0 <- results(dds0)[1:3,"log2FoldChange"]
  LFC <- results(dds1)[1:3,"log2FoldChange"]

  # filtered
  checkTrue(all(is.na(pval0)))
  # not filtered
  checkTrue(all(!is.na(pval[2:3])))

  # counts still the same
  checkTrue(all(counts(dds1)==counts(dds)))
  
  checkTrue(is.na(LFC[1]))
  # replaced, reduced LFC
  checkTrue(abs(LFC[2]) < abs(LFC0[2]))
  # replaced, LFC now zero
  checkTrue(LFC[3] == 0)
  
  idx <- which(!mcols(dds1)$replace)
  # the pvalue for those not replaced is equal
  checkEquals(results(dds1)$pvalue[idx], results(dds0)$pvalue[idx])

  # check that outlier filtering catches throughout range of mu
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100, m=12, interceptMean=seq(from=1,to=9,length=100), dispMeanRel=function(x) 4/x + .5)
  idx <- rep(rep(c(TRUE,FALSE),c(1,9)),10)
  counts(dds)[idx,1] <- 100000L
  dds <- DESeq(dds, minReplicatesForReplace=Inf)
  res <- results(dds)
  checkTrue(all(is.na(res$pvalue[idx])))
}
