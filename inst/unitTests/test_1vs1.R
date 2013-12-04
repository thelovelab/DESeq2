test_1vs1 <- function() {
  dds <- makeExampleDESeqDataSet(m=2)
  dds <- DESeq(dds)
  res <- results(dds)
}
