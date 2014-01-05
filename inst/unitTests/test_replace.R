test_replace <- function() {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100,m=12)
  counts(dds)[37,1] <- 100000L
  counts(dds)[54,7] <- 100000L
  dds1 <- DESeq(dds,minReplicatesForReplace=6)
  absLFC <- abs(results(dds1)[37,"log2FoldChange"])
  checkTrue(absLFC < .1)

  dds0 <- DESeq(dds,minReplicatesForReplace=7)
  idx <- which(!mcols(dds1)$replace)
  checkEquals(results(dds1)$pvalue[idx], results(dds0)$pvalue[idx])
}
