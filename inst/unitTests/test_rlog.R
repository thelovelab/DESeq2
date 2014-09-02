test_rlog <- function() {
  dds <- makeExampleDESeqDataSet(n=5000, m=20)
  rld <- rlog(dds, fast=TRUE)
}
