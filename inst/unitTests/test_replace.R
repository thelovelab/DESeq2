test_replace <- function() {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100, m=12)
  counts(dds)[1,] <- rep(0L, 12)
  counts(dds)[2,] <- c(100000L, rep(10L, 11))
  counts(dds)[3,] <- c(100000L, rep(0L, 11))
  dds1 <- DESeq(dds, minReplicatesForReplace=6)

  head(results(dds1),3)
  
  LFC <- results(dds1)[1:3,"log2FoldChange"]
  checkTrue(is.na(LFC[1]))
  checkTrue(abs(LFC[2]) < .2)
  checkTrue(LFC[3] == 0)
  
  dds0 <- DESeq(dds, minReplicatesForReplace=7)
  idx <- which(!mcols(dds1)$replace)
  checkEquals(results(dds1)$pvalue[idx], results(dds0)$pvalue[idx])
}
