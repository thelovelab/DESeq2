test_vst <- function() {
  dds <- makeExampleDESeqDataSet(n=100, m=4)
  design(dds) <- ~ 1
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, fitType="parametric")
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  dds <- estimateDispersions(dds, fitType="local")
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  dds <- estimateDispersions(dds, fitType="mean")
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)  
}
