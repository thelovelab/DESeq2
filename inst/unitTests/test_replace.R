test_replace <- function() {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100,m=12)
  counts(dds)[1,1] <- 100000L
  dds <- DESeq(dds,minReplicatesForReplace=6)
  checkTrue(abs(results(dds)[1,"log2FoldChange"]) < .1)
}
