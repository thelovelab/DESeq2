test_interactions <- function() {
  dds <- makeExampleDESeqDataSet()
  colData(dds)$group <- rep(letters[25:26],times=ncol(dds)/2)
  design(dds) <- ~ condition + condition:group
  dds <- DESeq(dds)
}
