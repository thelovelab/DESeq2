test_replace <- function() {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100, m=12)
  counts(dds)[37,] <- c(100000L, rep(10L, 11))
  counts(dds)[54,] <- c(100000L, rep(0L, 11))
  dds1 <- DESeq(dds, minReplicatesForReplace=6)
  absLFC <- abs(results(dds1)[37,"log2FoldChange"])
  checkTrue(absLFC < .2)

  dds0 <- DESeq(dds,minReplicatesForReplace=7)
  idx <- which(!mcols(dds1)$replace)
  checkEquals(results(dds1)$pvalue[idx], results(dds0)$pvalue[idx])
}
