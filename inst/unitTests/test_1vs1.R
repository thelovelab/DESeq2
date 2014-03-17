test_1vs1 <- function() {
  dds <- makeExampleDESeqDataSet(n=100, m=2)
  dds <- DESeq(dds)
  res <- results(dds)
}
